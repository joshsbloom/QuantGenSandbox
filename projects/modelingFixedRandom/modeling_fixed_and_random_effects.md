# Tutorial 2 - Heritability, Fixed and Random Effect Models

In our weekly journal club we come across an eyebrow-raising paper. In the paper, the authors claim to have found a novel anti-fungal compound. The authors show it inhibits the growth of many different kinds of fungi. However, the details about what exactly it is, whether it is fungal-specific, and how it works are frustratingly absent from the paper. The lab decides to request an aliquot of the compound, chemical X, from the authors to see if different strains of yeast are affected by it. If there are heritable differences in growth as a function of treatment with chemical X, perhaps a QTL mapping approach could help you find genetic variants that cause differential response to the compound. If you know the genes and the genetic variants in those genes that affect sensitivity or resistance to chemical X then maybe you can better understand what the compound is, how it is working and how resistance may evolve.

After a long (and ultimately unfocused and confusing) diatribe in lab meeting by Josh, where he seemed to go on endlessly about how much more could be achieved if we just worked together more effectively, the importance of randomization in experiment design, and on quantifying heritability before executing a genetic mapping experiment ... eight members of the lab come up with a plan to see if there are heritable differences in yeast growth in response to chemical X.  

Heriberto streaks out ten different yeast isolates and grows up overnight cultures of each and dilutes them back to 0.05 OD for innoculation. He also makes one working stock of YPD and YPD + Chemical X (YPD:X). Then the eight lab members each set up a subset of 48 wells across the plates, pipetting the appropriate media and innoculating it with the appropriate strain with their own pipettes. The four plates are grown overnight in the plate reader. Consistent with the disappointing nature of the current timeline, there was a power outage during the experiment, so all the lab is able to salvage is an endpoint OD reading for each plate.  

Josh says, *"No worries, there's plenty of information in the endpoint OD for you to estimate heritability from this! But why on earth were so many of you involved in this experiment setup, each using your own pipettes!?"* After the group explains the experiment design he says, *"Fine, whatever ... you kept track of everything relevant, time for some modeling!"*

---

## 1. Read in the data
- Read in the data from the GitHub repository [`joshsbloom/QuantGenSandbox`](https://github.com/joshsbloom/QuantGenSandbox):
- Read in the experiment design data [`projects/modelingFixedRandom/02_OD_design.csv`](/projects/modelingFixedRandom/02_OD_design.csv)
- Read in the plateR formatted endpoint OD data[`projects/modelingFixedRandom/02_OD_measurements.csv`](/projects/modelingFixedRandom/02_OD_measurements.csv).
- Use the **plateR** package to parse the plate data.
- Merge the measurements with the design into one `data.frame` for analysis.

---

## 2. Visualize the results

**a. Heatmaps**
- Use `ggplot2` to make **8×12 heatmaps** of OD for each plate.  
- Color wells using the **viridis** scale.  
- Are there any visually obvious plate- or position-dependent effects?  
- Make additional similarly structured heatmaps showing **experimenter by well position** and **treatment by well position**.

**b. Publication-ready 3-panel figure**
- Create sub-panels A–C, label them in the **top-left corner**. Use `patchwork` and `ggplot2`.
- Make **panel A twice the size** of B and C combined.
- Panel A: Violin plots of endpoint OD by **strain** (x-axis), separate violins for each treatment.  
  - Overlay raw data (jittered points).  
  - Show means for each group.  
  - Draw lines connecting means of YPD vs YPD:X within each strain.
- Panel B: Violin + jitter plots grouped by **experimenter**.  
- Panel C: Violin + jitter plots grouped by **plate**.  

**Questions:**  
- What do you conclude about strain, treatment, experimenter, and plate?  
- Does anything stand out visually?
- Why was it not ideal to have different people with different pipettes add the strain and media to each plate? What might you do differently next time? 

---

## 3. Fixed-effects linear model

- Fit a linear model:  
  ```
  OD ~ strain + treatment + plate + experimenter
  ```
- What is the **average effect of treatment**?  
- Use `Anova()` to test significance of each factor.  
- Next, add an **interaction** term:  
  ```
  OD ~ strain + treatment + strain:treatment + plate + experimenter
  ```
- Use a likelihood ratio test to decide if the inclusion of the interaction term is justified.
- Compare results to your visualization in 2b (panel A).

---

## 4. Broad-sense heritability in YPD:X

- Restrict analysis to YPD:X wells only.
- Estimate broad-sense heritability:
    ```
    H^2 = Var(Strain)/(Var(Strain)+Var(Residual))
    ```  
  - First **ignoring** experimenter and plate in the model.  
  - Then **including** experimenter and plate in the model.  
- How do the results differ?  
- Why does accounting for experimenter/plate change the estimates?  
- Why is heritability of this trait hard to interpret on its oown?

---

## 5. Mixed model with random effects
- Josh says plate and experimenter are better modeled as random effects here.
- Fit a mixed model:  
  ```
  OD ~ strain * treatment + (1|plate) + (1|experimenter)
  ```
- Compare this model (using **BIC**) to the fixed-effects model in step 3.  
- Was Josh right that plate and experimenter should be random effects? Beyond just differences (or no differences) in BIC, why or why not?

---

## 6. Normalization of OD in chemical X

- Get residuals from a model of OD with random effects of plate and experimenter (call this `r1`).  
- Compute mean residuals for each strain in **YPD** (call this `m1s`).  
- Subtract `m1s` from each strain’s residuals in **YPD:X**.  
- Retain only these corrected YPD:X values.  
- Make a violin plot of these normalized strain responses, including the jittered individually normalized data points in the plot.  

**Question:** What do you conclude about strain-level sensitivity to Chemical X after normalization?

---

## 7. Broad-sense heritability of the normalized trait

- Calculate H^2 for the normalized trait (as in step 6).  
- Reflect: **How can heritability be misleading?**  
**Question:** How would the information here and in step 6 guide a QTL mapping experiment you could perform.

---

