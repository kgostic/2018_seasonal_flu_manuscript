# Seasonal_flu

Created this repo to track R code for my seasonal flu project
Created after previous repo was not updating


# Branches
* Created branch back_in_time to change the way the likelihood function deals with years in which no cases are predicted.
* In master  the likelihood imposes an arbitrary penalty and likelihood is somewhat unstable
* In the new "back_in_time" strategy, go back to the closest year in which cases are predicted when no cases are predicted in the year of interest. Option to impose some small penalty as well.

