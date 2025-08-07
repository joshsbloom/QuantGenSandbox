# Yeast Growth in YNB + Glucose Media

## Background

You have two strains of yeast:  
- **BY**: Lab reference strain  
- **RM**: Vineyard isolate  

You want to determine whether there is a difference in how the yeast grow in **YNB + glucose media**.

---

## 1) Growth Curve Setup

- **120 µL cultures** of the strains were grown overnight to saturation in a 96-well plate called **`Plate080525_YNB_GLU`**.  
  - **BY**: Well `A01`  
  - **RM**: Well `A02`

- Use the **Biomek robot** to inoculate the two yeast strains into another 96-well plate called **`Plate080625perm_YNB_GLU`** to obtain replicate growth curves. 
- The **Biomek robot** needs a CSV file to tell it how to set up the 96-well plate; from which well to pipette from, which well to pipette into and how much volume.

- The 96-well plate format has rows `A`–`H` and columns `01`–`12`. **All numbers less than 10 need a preceding 0** (e.g., `B03`, `C09`) for the Biomek to process them correctly.

- **Goal**: Inoculate **1 µL** of saturated culture into **30 random wells** for each strain. The resulting plate will give you 30 replicate growth curves per strain.

- **Note**:  
  - **Edge wells** should contain **media blanks** to evaluate contamination and avoid unreliable curves due to evaporation.
  - If you pre-fill the plate manually, any wells **not listed** in the CSV will stay blank.

- **Biomek CSV format requirements**:
  - File must be **comma-delimited**
  - Use **DOS new line characters**
  - Include the following header:
    ```
    Source.Plate,Source.Wells,Strain,Target.Plate,Target.Wells,Volume
    ```

### Hints for Creating the File in R

- Use:
  - `data.frame()`
  - `sample()`
  - `set.seed(10)` specifying a specific seed will enable reproducible randomizations from the code
  - `paste0()`, `rep()`, `toupper()`, `sprintf()`, `sample()`
  - `write.table()` with:
    - `col.names`
    - `row.names`
    - `sep`
    - `quote`
    - `eol`

- **Extra credit**: Write the whole code in **4 lines**.

---

## 2) Analyzing the Growth Curves

You’ve completed the experiment. Is there a difference in how BY and RM grow?

### Step 1: Fit Growth Curves

- Read in the data from the GitHub repository [`joshsbloom/QuantGenSandbox`](https://github.com/joshsbloom/QuantGenSandbox):
  - plate reader growth data: [`projects/growthCurves/growthCurves_01.csv`](/projects/growthCurves/growthCurves_01.csv)
  - Biomek CSV file mapping strains to wells: [`projects/growthCurves/randomize_01.csv`](/projects/growthCurves/randomize_01.csv)

- Use:
  - `growth.gcFitSpline()` from the **QurvE** R package
  - `lubridate` to convert time into **decimal hours**

- Extract **doubling times** during log-phase growth.  
- Visualize a few of the growth curves and fits to confirm you understand the procedure and fitting looks reasonable.

---

### Step 2: Visualize Doubling Times Distribution by Strain

- Use `ggplot2` to make:
  - **Histogram** of doubling times by strain (with different colors)
  - **Violin plot** showing:
    - Individual points
    - Mean
    - 95% confidence interval of the mean

---

### Step 3: Statistical Test

- Use `t.test()` in R to test for differences in means.  
- **Question**: Should this be a **paired t-test**?

---

### Step 4: Permutation Test 

Josh is concerned about whether the assumptions of a t.test are valid here, so he recommends performing a **permutation test**.

- **Questions**: What assumptions does the permutation test make? What are the assumptions about a parametric t-test that Josh is concerned about?

- Generate a **null distribution** of t-statistics by **permuting strain labels** 1000 times and calculating the t-statistic each time
    - Visualize where the observed t-statistic is compared to this empirical null distribution
    - Calculate an **p-value** based on comparing the observed t-statistic for the non-permuted data vs the empirical null distribution

---

### Step 5: Power Calculations

Josh now asks:

> How small of a fitness difference (measured in doubling time) are we powered to detect between BY and an arbitrary new strain with a different doubling time?

Evaluate statistical power by simulation across this grid of parameters:

- **Replicates per strain**: 3, 6, 12, 24, 48, 96, 384  
- **Doubling time differences**: From 60 to 240 minutes, in **3-minute increments**  
- **Assumptions**:
  - BY mean doubling time = **90 min**
  - SD = **18 min**
  - Other strain has same SD

- For each combo of parameters:
  - Run 1000 simulations, simulating from two normal distributions with means being the expected doubling times, SD as specified above, and number of replicates per strain
  - Calculate power as the fraction of simulations out of 1000 per combo with t-test `p < 0.05`

- **Plot** using `ggplot2`:
  - **X-axis**: Difference in doubling time (minutes)  
  - **Y-axis**: Power  
  - **Line color/type**: Sample size (convert to factor for prettier colors)

- **Question**:  
  > If you want to detect a 10-minute difference in doubling time, how many replicates do you need?
