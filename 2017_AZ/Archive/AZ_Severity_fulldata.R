



Load the full data set
```{r}
raw.data = read.csv('AZ_2016_Hospitalization.csv')
raw.data$admitted = as.numeric(raw.data$day_admit != '0: no admit') # Indicator variable for admitted patients
```

Check whether the pattern holds for 1957 hospital admissions
Expect a lower proportion of H1N1 cases than H3N2 to occur in pre-1957 cohorts
```{r}
admit.table = wk.table*0; colnames(admit.table) = c('pre.1957', 'post.1957')
admit.table[1, 1] = nrow(subset(raw.data, admitted == 1 & MORBSPC == 'A/H1N1 (SWINE-ORIGIN)' & BYEAR < 1957))
admit.table[2, 1] = nrow(subset(raw.data, admitted == 1 & MORBSPC == 'A/H3' & BYEAR < 1957))
admit.table[1, 2] = nrow(subset(raw.data, admitted == 1 & MORBSPC == 'A/H1N1 (SWINE-ORIGIN)' & BYEAR >= 1957))
admit.table[2, 2] = nrow(subset(raw.data, admitted == 1 & MORBSPC == 'A/H3' & BYEAR >= 1957))
# Expect a higher fraction of H3N2 admits in the pre-1957 cohort than H1N1
prop.table(admit.table, 1)
c.test = chisq.test(admit.table)
c.test
c.test$expected
```


```{r echo = FALSE}
icu.table = wk.table*0; colnames(icu.table) = c('pre.1957', 'post.1957')
icu.table[1, 1] = nrow(subset(raw.data, ICU == 1 & MORBSPC == 'A/H1N1 (SWINE-ORIGIN)' & BYEAR < 1957))
icu.table[2, 1] = nrow(subset(raw.data, ICU == 1 & MORBSPC == 'A/H3' & BYEAR < 1957))
icu.table[1, 2] = nrow(subset(raw.data, ICU == 1 & MORBSPC == 'A/H1N1 (SWINE-ORIGIN)' & BYEAR >= 1957))
icu.table[2, 2] = nrow(subset(raw.data, ICU == 1 & MORBSPC == 'A/H3' & BYEAR >= 1957))
# Expect a higher fraction of H3N2 admits in the pre-1957 cohort than H1N1
icu.table
prop.table(icu.table, 1)
c.test = chisq.test(icu.table)
c.test
c.test$expected
icu.table-c.test$expected
````

Repeat for 1968 break
```{r echo = FALSE}
icu.table = wk.table*0; colnames(icu.table) = c('pre.1968', 'post.1968')
icu.table[1, 1] = nrow(subset(raw.data, ICU == 1 & MORBSPC == 'A/H1N1 (SWINE-ORIGIN)' & BYEAR < 1968))
icu.table[2, 1] = nrow(subset(raw.data, ICU == 1 & MORBSPC == 'A/H3' & BYEAR < 1968))
icu.table[1, 2] = nrow(subset(raw.data, ICU == 1 & MORBSPC == 'A/H1N1 (SWINE-ORIGIN)' & BYEAR >= 1968))
icu.table[2, 2] = nrow(subset(raw.data, ICU == 1 & MORBSPC == 'A/H3' & BYEAR >= 1968))
# Expect a higher fraction of H3N2 admits in the pre-1957 cohort than H1N1
icu.table
prop.table(icu.table, 1)
c.test = chisq.test(icu.table)
c.test
c.test$expected
icu.table-c.test$expected
```



Hospitalization rate?
```{r}
358+529 # H1N1 hosptialized
147+274 # H1N1 not hosptialized
(358+529)/(358+529+147+274) # H1N1 hosp rate

382+393 # H3N2 hospitalized
147+274 # H3N2 not hospitalized
(382+393)/(382+393+147+274) # H3N2 hosp rate
```