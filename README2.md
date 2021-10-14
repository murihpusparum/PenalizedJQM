## Using LQMM and penalized JQM to estimate Individual Reference Intervals

# Requirements

- To compile all the codes, R 4.0 or later is recommended.
- Install R package "quantreg", and "invgamma" (only for data generation used in the example) .

# Examples

- An example using generated data from Linear Mixed Models is available in "./examples" directory.
- A plot to illustrate the intervals on top of the data is produced. 
- The estimated intervals should be used to interpret the **future measurements**, not the existing data used in the estimation procedure.
- The empirical coverage is also computed i.e. the proportion of measurements inside the intervals. If we assume all the measurements are healthy (come from healthy individuals), we expect the empirical coverage to be close to the nominal level (e.g. 95%).
