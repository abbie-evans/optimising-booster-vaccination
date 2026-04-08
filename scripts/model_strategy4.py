import numpy as np
import pandas as pd
from scipy.stats import gamma, weibull_min
import random
from joblib import Parallel, delayed
from tqdm import tqdm
from scipy.integrate import quad

IC_df = pd.read_csv('data/Start_pop.csv',
                    skiprows=15, header=None, dtype=np.float64)
init = pd.read_csv('data/Start_pop.csv',
                   header=None, dtype=np.float64)
init.iloc[15] = init.iloc[15:].sum(axis=0)
init_cond = init.iloc[:16]
N = np.sum(init_cond, axis=1).tolist()
# Calculate the proportion of each age group in the population
prob = [i for i in N/np.sum(N)]

frac_pop_over75 = ((1/np.sum(np.asarray(IC_df))) * np.sum(np.asarray(IC_df),axis=1))
hosp_risk = [0.01125507, 0.01092597, 0.00612118, 0.00493643, 0.00394915, 0.0027644,
            0.0041466, 0.00592372, 0.00789829, 0.0106627, 0.01072852, 0.01020196,
            0.0137562, 0.01625732, 0.01599405, 0.01711297, 0.01862681, 0.01737625,
            0.0151384, 0.01276891, 0.01276891]
death_risk = [0.00148776, 0.00104143, 0.01413369, 0.00773634, 0.00892654, 0.01874574,
            0.01740676, 0.01874574, 0.02796984, 0.03109413, 0.04745945, 0.08450461,
            0.14565143, 0.13702244, 0.24592627, 0.33712578, 0.38964361, 0.51699563,
            0.86007244, 1, 1]
beta = pd.read_csv('data/UK_Risks.csv')['susceptibility'].tolist()
# -------------------------------------------------------------
# Calculate gamma distributions for hospitalisation and death
mean_h = 0.00001
mean_d = 10
std_dev = 12.1
shape_h = (mean_h/std_dev)**2
scale_h = std_dev**2/mean_h
shape_d = (mean_d/std_dev)**2
scale_d = std_dev**2/mean_d
deltaIH = gamma.pdf(np.arange(1, 31), a=shape_h, scale=scale_h)
deltaIH /= np.sum(deltaIH)
deltaHD = gamma.pdf(np.arange(1, 31), a=shape_d, scale=scale_d)
deltaHD /= np.sum(deltaHD)
# -------------------------------------------------------------
M = pd.read_csv('data/UK_Contacts.csv', header=None)
M_list = [M.iloc[i,:].tolist() for i in range(16)]
d = pd.read_csv('data/UK_Risks.csv')['symptom_risk'].tolist()
d = np.array(d[:15] + [np.sum(np.multiply(d[15:], frac_pop_over75))])

def l_k(k, shape, scale):
    integrand = lambda u: (1 - abs(u - k)) * gamma.pdf(u, a=shape, scale=scale)
    result, _ = quad(integrand, k - 1, k + 1)
    return result

# Latent period distribution
l_values = [l_k(k, shape=3, scale=5/3) for k in range(2, 40)]
# l_1 is 1-sum of l_values
l_1 = 1 - sum(l_values)
l_period = [l_1] + l_values

# Infectious period distribution
l_values = [l_k(k, shape=3, scale=9/3) for k in range(2, 40)]
# l_1 is 1-sum of l_values
l_1 = 1 - sum(l_values)
inf_period = [l_1] + l_values

# Delay to hospital distribution
def h_k(k, shape, scale):
    integrand = lambda u: (1 - abs(u - k)) * weibull_min.pdf(u, c=shape, scale=scale)
    result, _ = quad(integrand, k - 1, k + 1)
    return result
h_values = [h_k(k, shape=1.4, scale=8.4) for k in range(2, 31)]
h_1 = 1 - sum(h_values)
deltaIH = [h_1] + h_values

# Delay to death distribution
d_values = [l_k(k, shape=shape_d, scale=scale_d) for k in range(2, 31)]
d_1 = 1 - sum(d_values)
deltaHD = [d_1] + d_values

parameters = {
    # Disease transmission and progression rates
    "beta": np.array(beta[:15] + [np.sum(np.multiply(beta[15:], frac_pop_over75))]), # Infection rate parameter for each age group
    "omega": 0.304,
    "prob_symptoms_by_age": d,
    "d": d, # Probability of symptomatic infection by age group

    # Hospitalisation parameters
    "pIH": np.array(hosp_risk[:15] + [np.sum(np.multiply(hosp_risk[15:], frac_pop_over75))]), # Probability of infection leading to hospitalisation by age group
    "deltaIH": deltaIH,
    'pHD': np.array(death_risk[:15] + [np.sum(np.multiply(death_risk[15:], frac_pop_over75))]), # Probability of hospitalisation leading to death by age group
    'deltaHD': deltaHD,

    # Boosting parameters
    "num_boost": 2000, # Number of individuals to boost daily

    # Contact matrix for force of infection
    "M_list": M_list,
}

def protection_from_infection(multiplier):
    """
    Calculate the protection from infection curves.
    """
    n_i0 = 1.13
    n_v0 = 1.13
    pi1 = -np.log(2)/35
    pi2 = -np.log(2)/1000
    ts = 75
    k = 2.5
    n50 = 0.091
    tau = np.linspace(0, 2000, 2000)

    y = lambda t: n_v0*((np.exp(pi1*t+pi2*ts)+np.exp(pi2*t+pi1*ts))/(np.exp(pi1*ts)+np.exp(pi2*ts)))
    i = lambda t: n_i0*(np.exp(pi1*t+pi2*ts)+np.exp(pi2*t+pi1*ts)/(np.exp(pi1*ts)+np.exp(pi2*ts)))
    epsilon = lambda n: 1/(1+np.exp(-k*(np.log10(n/5.1)-np.log10(n50)))) # Existing booster
    epsilon2 = lambda n: 1/(1+np.exp(-k*(np.log10(n*multiplier/5.1)-np.log10(n50)))) # Variant-adapted booster
    epsilon3 = lambda n: 1/(1+np.exp(-k*(np.log10(n*3/5.1)-np.log10(n50)))) # Variant infection

    return epsilon(y(tau)), epsilon2(y(tau)), epsilon3(i(tau))

def calc_risk_of_hospitalisation(multiplier):
    n_i0 = 1.13
    n_v0 = 1.13
    pi1 = -np.log(2)/35
    pi2 = -np.log(2)/1000
    ts = 75
    k = 2.5
    n50 = 0.021
    tau = np.linspace(0, 2000, 2000)

    y = lambda t: n_v0*((np.exp(pi1*t+pi2*ts)+np.exp(pi2*t+pi1*ts))/(np.exp(pi1*ts)+np.exp(pi2*ts)))
    i = lambda t: n_i0*(np.exp(pi1*t+pi2*ts)+np.exp(pi2*t+pi1*ts)/(np.exp(pi1*ts)+np.exp(pi2*ts)))
    epsilon = lambda n: 1/(1+np.exp(-k*(np.log10(n/5.1)-np.log10(n50)))) # Existing booster
    epsilon2 = lambda n: 1/(1+np.exp(-k*(np.log10(n*multiplier/5.1)-np.log10(n50)))) # Variant-adapted booster
    epsilon3 = lambda n: 1/(1+np.exp(-k*(np.log10(n*3/5.1)-np.log10(n50)))) # Variant infection

    return epsilon(y(tau)), epsilon2(y(tau)), epsilon3(i(tau))

def calculate_susc(t, boosted, boost_time, recovered_time, entry_time, susceptibility, hosp_protection, vaccine):
    # Calculate susceptibility based on the population's vaccination status and time since last boost
    
    not_boosted_not_recovered = (~boosted) & (np.isnan(recovered_time))
    susceptibility[not_boosted_not_recovered] = 1 - protect_infection[0][np.abs(t + entry_time[not_boosted_not_recovered]).astype(int)]
    hosp_protection[not_boosted_not_recovered] = risk_of_hospitalisation[0][np.abs(t + entry_time[not_boosted_not_recovered]).astype(int)]

    # Condition 2: Boosted, not recovered
    boosted_not_recovered = (boosted) & (np.isnan(recovered_time)) & (vaccine==1)
    boost_time_diff = (t - boost_time[boosted_not_recovered]).astype(int)
    susceptibility[boosted_not_recovered] = 1 - protect_infection[1][boost_time_diff]
    hosp_protection[boosted_not_recovered] = risk_of_hospitalisation[1][boost_time_diff]

    # Condition 2b: Boosted, not recovered, existing vaccine
    boosted_not_recovered_existing = (boosted) & (np.isnan(recovered_time)) & (vaccine==0)
    boost_time_diff = (t - boost_time[boosted_not_recovered_existing]).astype(int)
    susceptibility[boosted_not_recovered_existing] = 1 - protect_infection[0][boost_time_diff]
    hosp_protection[boosted_not_recovered_existing] = risk_of_hospitalisation[0][boost_time_diff]

    # Condition 3: Not boosted, recovered
    not_boosted_recovered = (~boosted) & (~np.isnan(recovered_time))
    recovered_time_diff = (t - recovered_time[not_boosted_recovered]).astype(int)
    susceptibility[not_boosted_recovered] = 1 - protect_infection[2][recovered_time_diff]
    hosp_protection[not_boosted_recovered] = risk_of_hospitalisation[2][recovered_time_diff]

    # Condition 4: Boosted and recovered
    boosted_recovered = (boosted) & (~np.isnan(recovered_time)) & (vaccine==1)
    boost_time_diff = (t - boost_time[boosted_recovered]).astype(int)
    recovered_time_diff = (t - recovered_time[boosted_recovered]).astype(int)
    options = [protect_infection[1][boost_time_diff], protect_infection[2][recovered_time_diff]]
    max_option = [0 if pair[0] >= pair[1] else 1 for pair in zip(*options)]
    susceptibility[boosted_recovered] = 1 - np.array([options[max_option[i]][i] for i in range(len(max_option))])
    options = [risk_of_hospitalisation[1][boost_time_diff], risk_of_hospitalisation[2][recovered_time_diff]]
    max_option = [0 if pair[0] >= pair[1] else 1 for pair in zip(*options)]
    hosp_protection[boosted_recovered] = np.array([options[max_option[i]][i] for i in range(len(max_option))])

    # Condition 4b: Boosted and recovered with existing vaccine
    boosted_recovered = (boosted) & (~np.isnan(recovered_time)) & (vaccine==0)
    boost_time_diff = (t - boost_time[boosted_recovered]).astype(int)
    recovered_time_diff = (t - recovered_time[boosted_recovered]).astype(int)
    options = [protect_infection[0][boost_time_diff], protect_infection[2][recovered_time_diff]]
    max_option = [0 if pair[0] >= pair[1] else 1 for pair in zip(*options)]
    susceptibility[boosted_recovered] = 1 - np.array([options[max_option[i]][i] for i in range(len(max_option))])
    options = [risk_of_hospitalisation[0][boost_time_diff], risk_of_hospitalisation[2][recovered_time_diff]]
    max_option = [0 if pair[0] >= pair[1] else 1 for pair in zip(*options)]
    hosp_protection[boosted_recovered] = np.array([options[max_option[i]][i] for i in range(len(max_option))])

def calculate_force_of_infection(init_susc, age_group, status):
    omega = parameters['omega']
    M_list = np.array(parameters['M_list']).T
    tau = 0.255
    d = parameters['d']
    R = np.zeros((16, 16))
    for i in range(16):
        for j in range(16):
            R[i, j] = (parameters['beta'][j]*parameters['omega']*M_list[i, j]*init_susc[i]*9)*(tau*(1-d[j]) + d[j])
    eigenvalues = np.linalg.eigvals(R)
    R0 = np.real(max(eigenvalues))
    beta_new = 1.5/R0*parameters['beta'] # Set Re
    N = np.bincount(age_group, minlength=16)
    infected_age_groups = age_group[status == 2]
    asymptomatic_age_groups = age_group[status == 3]
    infected_counts = np.bincount(infected_age_groups, minlength=16)
    asymptomatic_counts = np.bincount(asymptomatic_age_groups, minlength=16)
    total_infected = infected_counts + tau * asymptomatic_counts
    force_of_infection = beta_new * omega * (M_list @ total_infected) / N
    return force_of_infection

def boost_people(current_time, status, boosted, boost_time, to_boost_time, dead, pop, vaccine):
    to_boost = ((status[pop] == 0) | (status[pop] == 1)) & (~boosted[pop]) & (~dead[pop]) & (to_boost_time[pop] <= current_time)
    # Select individuals to boost based on the number of boosts available
    num_boost = parameters["num_boost"]
    # if we can vaccination with vaccine=0 and vaccine=1 at the same time (old and new), prioritise the new vaccine (vaccine=1)
    # if in to_boost list, there are individuals with vaccine=1, boost them first
    to_boost_indices = pop[to_boost]
    if vaccine[to_boost_indices].sum() > 0:
        # reorder to_boost so that individuals with vaccine=1 are at the front
        vaccine_1_indices = to_boost_indices[vaccine[to_boost_indices] == 1]
        vaccine_0_indices = to_boost_indices[vaccine[to_boost_indices] == 0]
        reordered_indices = np.concatenate((vaccine_1_indices, vaccine_0_indices))
        to_boost_indices = reordered_indices
    else:
    # Find indices of pop where to_boost is True
        to_boost_indices = pop[to_boost]
    # Select up to num_boost individuals from to_boost_indices
    selected_indices = to_boost_indices[:num_boost]
    # print(vaccine[selected_indices], "people boosted at time", current_time)
    boost_time[selected_indices] = current_time
    boosted[selected_indices] = True

def filter_population(age_group):
    # sort the population by increasing age group and choose 80% of each age group
    age_groups = np.arange(16)
    age_group_counts = np.bincount(age_group, minlength=16)
    filtered_population = []
    over_50 = 0
    for age in age_groups:
        count = int(age_group_counts[age] * 0.8)
        over_50 += count if age >= 10 else 0
        indices = np.where(age_group == age)[0]
        if len(indices) > count:
            selected_indices = np.random.choice(indices, size=count, replace=False)
        else:
            selected_indices = indices
        filtered_population.extend(selected_indices)
    filtered_population = np.array(filtered_population)
    return filtered_population, over_50

def run(delay, old_start, sim_num=0):
    base_seed = 42
    random.seed(base_seed)
    np.random.seed(base_seed)
    ##########################
    # Initialise population
    ##########################
    status = np.zeros(100000)
    age_group = np.random.choice(np.arange(0, 16), size=100000, p=prob)
    entry_time = np.random.randint(0, 711, size=100000) # Random entry time for each individual
    latent_period = np.full(100000, np.nan)
    infectious_period = np.full(100000, np.nan)
    hospital_entry = np.full(100000, np.nan)
    death_entry = np.full(100000, np.nan)
    dead = np.zeros(100000, dtype=bool)
    boosted = np.zeros(100000, dtype=bool)
    boost_time = np.full(100000, np.nan)
    to_boost_time = np.full(100000, delay)
    recovered_time = np.full(100000, np.nan)
    susceptibility = np.ones(100000)
    hosp_protection = np.zeros(100000)
    exposed_entry_time = np.full(100000, np.nan)
    sympt_entry_time = np.full(100000, np.nan)
    hosp_prob = parameters['pIH'][age_group]
    vaccine = np.ones(100000, dtype=int)

    exposed_indices = np.random.choice(np.arange(100000), size=100, replace=False)
    status[exposed_indices] = 1 # STATUS_EXPOSED
    latent_period[exposed_indices] = np.random.choice(range(1, 40), size=100, p=l_period)
    exposed_entry_time[exposed_indices] = 0
    calculate_susc(0, boosted, boost_time, recovered_time, entry_time, susceptibility, hosp_protection, vaccine)
    hosp_prob[exposed_indices] = parameters['pIH'][age_group[exposed_indices]] * (1 - hosp_protection[exposed_indices]) / susceptibility[exposed_indices]

    base_seed = 42 + sim_num
    random.seed(base_seed)
    np.random.seed(base_seed)

    new_death = []
    age_death = []
    new_infect = []
    pop, over_50 = filter_population(age_group)
    vaccine[pop[:over_50]] = 0 # Set the first x people in the population to get the existing vaccine
    to_boost_time[pop[:over_50]] = old_start # Set the boost time for the first x people to the old start time

    # remaining population to get the new vaccine
    to_boost_time[pop[over_50:]] = delay
    vaccine[pop[over_50:]] = 1

    for t in range(365):
        if t == 0:
            calculate_susc(t, boosted, boost_time, recovered_time, entry_time, susceptibility, hosp_protection, vaccine) # Initialise susceptibility
            init_susc = np.zeros(16)
            for age in range(16):
                init_susc[age] = np.mean(susceptibility[age_group == age])

        force_of_infection = calculate_force_of_infection(init_susc, age_group, status)
        if t> 0:
            calculate_susc(t, boosted, boost_time, recovered_time, entry_time, susceptibility, hosp_protection, vaccine)

        susceptible_indices = np.where(status == 0)[0]
        exposed_indices = np.where(status == 1)[0]
        symptomatic_indices = np.where(status == 2)[0]
        asymptomatic_indices = np.where(status == 3)[0]

        # Susceptible -> Exposed
        prob_infection = 1 - np.exp(-susceptibility[susceptible_indices]*force_of_infection[age_group[susceptible_indices]])
        susceptible_to_exposed = (np.random.rand(len(susceptible_indices)) < prob_infection) & (np.isnan(death_entry[susceptible_indices]))
        new_infect.append(susceptible_to_exposed.sum())
        status[susceptible_indices[susceptible_to_exposed]] = 1 # STATUS_EXPOSED
        latent_period[susceptible_indices[susceptible_to_exposed]] = np.random.choice(range(1, 40), size=np.sum(susceptible_to_exposed), p=l_period)
        exposed_entry_time[susceptible_indices[susceptible_to_exposed]] = t

        # for the susceptible_to_exposed individuals, calculate the hospital entry
        pIH = (1-hosp_protection[susceptible_indices[susceptible_to_exposed]]) / susceptibility[susceptible_indices[susceptible_to_exposed]]
        hosp_prob[susceptible_indices[susceptible_to_exposed]] = pIH * parameters['pIH'][age_group[susceptible_indices[susceptible_to_exposed]]]

        # Exposed -> Symptomatic
        prob_symptoms = parameters['prob_symptoms_by_age'][age_group[exposed_indices]]
        exposed_to_symptomatic = (np.random.rand(len(exposed_indices)) < prob_symptoms) & (t-exposed_entry_time[exposed_indices] > latent_period[exposed_indices])
        exposed_to_asymptomatic = (~exposed_to_symptomatic) & (t-exposed_entry_time[exposed_indices] > latent_period[exposed_indices])
        status[exposed_indices[exposed_to_symptomatic]] = 2 # STATUS_SYMPTOMATIC
        status[exposed_indices[exposed_to_asymptomatic]] = 3 # STATUS_ASYMPTOMATIC
        infectious_period[exposed_indices[exposed_to_symptomatic]] = np.random.choice(range(1, 40), size=np.sum(exposed_to_symptomatic), p=inf_period)
        infectious_period[exposed_indices[exposed_to_asymptomatic]] = np.random.choice(range(1, 40), size=np.sum(exposed_to_asymptomatic), p=inf_period)
        sympt_entry_time[exposed_indices[exposed_to_symptomatic]] = t
        sympt_entry_time[exposed_indices[exposed_to_asymptomatic]] = t

        hosp_entry = (np.random.rand(len(exposed_indices[exposed_to_symptomatic])) < hosp_prob[exposed_indices[exposed_to_symptomatic]]) & (np.isnan(hospital_entry[exposed_indices[exposed_to_symptomatic]]))
        hospital_entry[exposed_indices[exposed_to_symptomatic][hosp_entry]] = t + np.random.choice(np.arange(1, 31), size=np.sum(hosp_entry), p=deltaIH) +1
        # for those in hospital, calculate the death entry
        death_prob = parameters['pHD'][age_group[exposed_indices[exposed_to_symptomatic][hosp_entry]]]
        death_entry_ = (np.random.rand(len(exposed_indices[exposed_to_symptomatic][hosp_entry])) < death_prob) & (np.isnan(death_entry[exposed_indices[exposed_to_symptomatic][hosp_entry]]))
        death_entry[exposed_indices[exposed_to_symptomatic][hosp_entry][death_entry_]] = hospital_entry[exposed_indices[exposed_to_symptomatic][hosp_entry][death_entry_]] + np.random.choice(np.arange(1, 31), size=np.sum(death_entry_), p=deltaHD) +1

        # Symptomatic -> Susceptible
        symptomatic_to_susceptible = (t - sympt_entry_time[symptomatic_indices] > infectious_period[symptomatic_indices])
        asymptomatic_to_susceptible = (t - sympt_entry_time[asymptomatic_indices] > infectious_period[asymptomatic_indices])
        status[symptomatic_indices[symptomatic_to_susceptible]] = 0 # STATUS_SUSCEPTIBLE
        status[asymptomatic_indices[asymptomatic_to_susceptible]] = 0 # STATUS_SUSCEPTIBLE
        recovered_time[symptomatic_indices[symptomatic_to_susceptible]] = t
        recovered_time[asymptomatic_indices[asymptomatic_to_susceptible]] = t
        entry_time[symptomatic_indices[symptomatic_to_susceptible]] = t
        entry_time[asymptomatic_indices[asymptomatic_to_susceptible]] = t
        
        death = (death_entry == t)
        # count number of deaths
        num_deaths = death.sum()
        # count number of deaths in each age group
        age_groups = np.arange(16)
        age_group_deaths = np.zeros(16, dtype=int)
        for age in age_groups:
            age_group_deaths[age] = np.sum(death[age_group == age])
        new_death.append(num_deaths)
        dead[death] = True
        age_death.append(age_group_deaths)
        status[death] = 0 # Move out of infectious status

        boost_people(t, status, boosted, boost_time, to_boost_time, dead, pop, vaccine)

    return new_infect, new_death, age_death

protect_infection = protection_from_infection(multiplier=2) # Set value of sigma
risk_of_hospitalisation = calc_risk_of_hospitalisation(multiplier=2) # Set value of sigma

for delay in np.linspace(0, 180, 181, dtype=int):
    inf = []
    death = []
    avg_age_deaths = []
    for old_start in range(0, delay + 1):
    # Run the simulation with the specified delay
        all_results = Parallel(n_jobs=-1)(delayed(run)(delay, old_start, _) for _ in tqdm(range(100)))
        all_infect, all_deaths, all_age_deaths = zip(*all_results)
        death.append(np.mean(all_deaths, axis=0))
        age_deaths = np.mean(all_age_deaths, axis=0)
        avg_age_deaths.append(age_deaths.sum(axis=0))
        inf.append(np.mean(all_infect, axis=0))

# Save the results to a CSV file
    df = pd.DataFrame(inf)
    df.to_csv(f'infections_s4_{delay}.csv', header=None, index=False)
    df = pd.DataFrame(death)
    df.to_csv(f'deaths_s4_{delay}.csv', header=None, index=False)
    age_deaths_df = pd.DataFrame(avg_age_deaths)
    age_deaths_df.to_csv(f'age_deaths_s4_{delay}.csv', header=None, index=False)
