#### PROJECT: Genomic offsets and demographic trajectories of Mimulus cardinalis populations during extreme drought
#### PURPOSE OF THIS SCRIPT: Test whether the dnorm estimation of growth is sufficient
#### AUTHOR: Amy Angert
#### DATE LAST MODIFIED: 20230224

## Note: Line 80 of script 02_lambda_IPM.R had a warning note that growth standard deviation method might be insufficient. This refers to an email conversation between Seema and Megan Demarche dating back to March 2019, in which Megan explained how using a cdf function could be better if variance in growth is small &/or if using a small number of mesh points. Megan's recommendation as of Feb 16, 2023 was to test if lambda is stabilizing as a function of n mesh points, and if so, then dnorm is ok to use. If not, then switch to cdf.

## Note: Outside of this script, I've manually set mesh points on line 44 of integral_projection_model.R to n=20 up to n=100 (our default) and rerun lines 217-365 of script 02_lambda_IPM.R for each n to calculate lambda at different mesh points. Here we will join those output files and examine lambda ~ n mesh.

colnames(site.info.20)[9] = "lambda.20"
colnames(site.info.30)[9] = "lambda.30"
colnames(site.info.40)[9] = "lambda.40"
colnames(site.info.50)[9] = "lambda.50"
colnames(site.info.60)[9] = "lambda.60"
colnames(site.info.70)[9] = "lambda.70"
colnames(site.info.80)[9] = "lambda.80"
colnames(site.info.90)[9] = "lambda.90"
colnames(site.info.100)[9] = "lambda.100"

lambdas.mesh <- left_join(left_join(left_join(left_join(left_join(left_join(left_join(left_join(site.info.20, site.info.30), site.info.40), site.info.50), site.info.60), site.info.70), site.info.80), site.info.90), site.info.100)

lambdas.long <- lambdas.mesh %>% 
  pivot_longer(cols=starts_with("lambda"), 
               names_to="mesh", 
               values_to="lambda") %>% 
  mutate(mesh.num = as.numeric(str_remove(mesh, "lambda.")))

ggplot(data=filter(lambdas.long, lambda>0), aes(x=mesh.num, y=lambda, group=Site, col=Site)) +
  geom_point() + 
  geom_smooth(fill=NA)
# Note: I don't know why the fitted lines sometimes have too low of a y-intercept, but regardless it is clear that lambda is stable across mesh points for all sites and we do not need to switch to cdf estimation

