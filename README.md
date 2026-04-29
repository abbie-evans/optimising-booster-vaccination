Data and Python code accompanying Evans et al., "Optimising booster vaccination against SARS-CoV-2 variants: When should a variant-adapted vaccine be developed?"

## Main scripts

### Data generation
- `model.py` \
  Simulates booster vaccination with either the deployment of an existing vaccine (with `multiplier=1`) available from the beginning of the simulation or the deployment of a variant-adapted vaccine (with `multiplier>1`) when it becomes available, where `multiplier` represents the relative vaccine efficacy within the code. Both of these vaccines are deployed in decreasing age order (starting with individuals aged 75+, then 70-74, and so on).
- `model_init_infect_50` and `model_init_infect_200` simulate booster vaccination as above, but with lower or higher numbers of initially infected individuals.
- `model_strategy3.py` \
  Simulates booster vaccination with the existing vaccine deploying $y$ doses. When the variant-adapted vaccine becomes available, vaccinate the remaining 80,000 - $y$ of the vaccine-willing population (based on the assumed vaccine uptake of 80%).

- `model_strategy4.py` \
  Simulates booster vaccination with the existing vaccine deploying $y$ doses in increasing age order (starting with individuals aged 0-4, then 5-9, and so on). When the variant-adapted vaccine becomes available, vaccinate the remaining 80,000 - $y$ of the vaccine-willing population (based on the assumed vaccine uptake of 80%).

### Plotting

To reproduce the figures, run the following scripts:

- Figure 1:
- Figure 2:
- Figure 3: `vaccine_efficacy_delay_heatmap.py`
- Figure 4: `alt_strategies.py`
- Figure S1:
- Figure S2: `diff_init_infect_heatmap.py`
- Figure S3: `vaccine_efficacy_delay_heatmap.py`
- Figure S4: `alt_strategies.py`
- Figure S5:
- Figure S6: `contact_matrix.py`

Note that all dependencies are specified in the `requirements.txt` file.

Reproduced figures will be saved into the `Figures` directory.
