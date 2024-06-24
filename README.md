# FAVAR Implementation in R

Welcome to the FAVAR (Factor Augmented Vector Auto-Regression) R implementation repository! This repository provides a comprehensive implementation of a FAVAR model for forecast prediction tasks in R.

## Table of Contents

- Introduction
- Installation
- Usage
  - Function Parameters
  - Example
- Dependencies
- Contributing
- License

## Introduction

FAVAR (Factor Augmented Vector Auto-Regression) is an advanced econometric model that extends the standard VAR (Vector Auto-Regression) model by incorporating factors extracted from a large dataset. This approach can improve forecast accuracy by leveraging a broader information set.

This implementation allows you to fit a FAVAR model and make forecasts using your own dataset in a seamless manner. 

## Installation

To get started, you need to clone the repository and install the required R packages:

```{r, engine='bash', count_lines}
git clone https://github.com/PetrGarm/FAVAR.git
cd FAVAR/Code/my experiments/FAVAR_CV.R
```

Next, ensure you have the necessary R packages installed:

```{r, engine='R', count_lines}
install.packages(c("tidyverse", "forecast", "vars", "lmtest"))
````

The main function provided by this repository is fore_FAVAR, which takes in your dataset and parameters and returns the forecast results.

### Function Parameters

- `X`: A matrix or data frame containing the predictors.
- `Y`: A matrix or data frame containing the dependent variable(s).
- `K`: The number of factors to extract from X.
- `y_name`: The name of the dependent variable in Y to be forecasted.
- `h`: The forecast horizon.
- `y`: (Optional and must not be changed probably) A time series object of the dependent variable. Defaults to Y[,y_name].
- `use_VAR`: (Optional) Logical flag to use standard VAR without factors if TRUE. Defaults to FALSE.

### Example

```{r, engine='R', count_lines}
# Sample usage of the fore_FAVAR function
source("FAVAR_CV.R")


# Define your data
X <- ...  # Your predictors matrix
Y <- ...  # Your dependent variables matrix
K <- 3    # Number of factors to extract
y_name <- "target_variable"  # The name of the dependent variable to forecast
h <- 12   # Forecast horizon

# Run FAVAR
forecast_result <- fore_FAVAR(X, Y, K, y_name, h)

# Display forecast results
print(forecast_result)

```

The following R packages are required to run the FAVAR implementation:

- `tidyverse`
- `forecast`
- `vars`
- `lmtest`

Ensure they are installed and loaded into your R environment before running the fore_FAVAR function.

## Contributing

Contributions to this repository are welcome. If you have any improvements, bug fixes, or new features, feel free to open a pull request or issue.

## License

This project is licensed under the MIT License. See the LICENSE file for more details.
...
