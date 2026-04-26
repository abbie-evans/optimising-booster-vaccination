Data and Python code accompanying Evans et al., "Optimising booster vaccination against SARS-CoV-2 variants: When should a variant-adapted vaccine be developed?"

## Main scripts

### Data generation
- `model.py` \
  Simulates booster vaccination with either the deployment of an existing vaccine (with `multiplier=1`) available from the beginning of the simulation or the deployment of a variant-adapted vaccine (with `multiplier>1`) when it becomes available, where `multiplier` represents the relative vaccine efficacy within the code. Both of these vaccines are deployed in decreasing age order (starting with individuals aged 75+, then 70-74, and so on).
- `model_strategy3.py` \
  Simulates booster vaccination with the existing vaccine deploying $y$ doses. When the variant-adapted vaccine becomes available, vaccinate the remaining 80,000 - $y$ of the vaccine-willing population (based on the assumed vaccine uptake of 80%).

- `model_strategy4.py` \
  Simulates booster vaccination with the existing vaccine deploying $y$ doses in increasing age order (starting with individuals aged 0-4, then 5-9, and so on). When the variant-adapted vaccine becomes available, vaccinate the remaining 80,000 - $y$ of the vaccine-willing population (based on the assumed vaccine uptake of 80%).

### Plotting

To reproduce the figures, run the following scripts:

- Figure 1:
- Figure 2:
- Figure 3:
- Figure 4:
- Figure S1:
- Figure S2:
- Figure S3:
- Figure S4:
- Figure S5:
- Figure S6:

Note that all dependencies are specified in the `requirements.txt` file.

Reproduced figures will be saved into the `Figures` directory.
