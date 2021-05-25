# hw2
There are three function in this hw2 package:

olve_ols(): a function solves a linear system using Gauss-Seidel or Jacobi method,
allows user to specify how many cores to use for parallel implementation.

algo_leverage(): a function implements algorithmic leveraging for linear 
regression using uniform and leverage score based subsampling of rows.

elnet_coord(): a function fits elastic net to data using coordinate descent algorithm.

## Download and install the package

```{r}
library(devtools)
install_github("https://github.com/AnnyShen55/hw2", force = TRUE)
library(hw2)
```

Then you will be able to use these three functions in this package.



