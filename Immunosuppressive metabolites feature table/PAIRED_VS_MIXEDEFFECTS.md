# Paired Analysis vs Mixed-Effects Models: A Practical Guide

## TL;DR - Quick Decision

| **Scenario** | **Recommendation** |
|--------------|-------------------|
| Exactly 2 timepoints per patient | **Use 06** (paired t-test/Wilcoxon) |
| Exactly 3 timepoints per patient | **Use 06b** (Friedman + pairwise) |
| Variable or >3 timepoints per patient | **Use 07** (mixed-effects) |
| Very small sample (n<10 patients) | **Use 06/06b** (simpler, more robust) |
| Moderate to large sample | **Use both 06 and 07** (complementary) |
| Want to test time×rejection interaction | **Use 07** (mixed-effects) |
| Want individual recovery patterns | **Use 06b** (reversibility classification) |

---

## What's the Difference?

### Paired Analysis (Scripts 06 and 06b)

**Philosophy:** Compare specific timepoint pairs within patients

**Approach:**
- Select specific timepoints (e.g., pre-2R vs active-2R)
- Require matched pairs or triplets
- Use simple paired tests (Wilcoxon, paired t-test, Friedman)
- Intuitive interpretation: "what changed from A to B?"

**Example:**
```
Patient 1: pre=10, active=15  → difference = +5
Patient 2: pre=12, active=18  → difference = +6
Patient 3: pre=8,  active=14  → difference = +6
Mean change: +5.67, p = 0.01 (paired Wilcoxon)
```

### Mixed-Effects Models (Script 07)

**Philosophy:** Model overall trajectory accounting for patient-specific baselines

**Approach:**
- Use ALL available timepoints from ALL patients
- Fit a statistical model with fixed and random effects
- Account for within-patient correlation
- Test complex hypotheses (time, ACR, ACR×time)

**Example:**
```
Model: metabolite ~ ACR_status + Time + (1|Patient)
         - Fixed effects: ACR and Time (population-level)
         - Random effects: Each patient has their own baseline
         
Result: ACR effect = +6.2 (SE=1.5), p = 0.001
        Accounts for fact that some patients naturally have higher/lower levels
```

---

## Key Conceptual Differences

### 1. Data Requirements

**Paired Analysis:**
- ✅ Works with unbalanced data BUT requires specific timepoints
- ❌ Patients without BOTH/ALL timepoints are excluded
- Example: For 06, need both pre AND active; missing either = excluded

**Mixed-Effects:**
- ✅ Uses ALL available data from ALL patients
- ✅ Patients with any number of timepoints can contribute
- ✅ Handles missing data naturally (maximum likelihood)
- Example: Patient with 1 sample still contributes to estimating population effects

### 2. Statistical Power

**Paired Analysis:**
```
Power depends on: Number of complete pairs
If n=20 patients, but only 10 have both pre+active → n=10 for analysis
```

**Mixed-Effects:**
```
Power depends on: Total number of observations + number of patients
If n=20 patients with 3 samples each → n=60 observations
Even patients with incomplete data contribute
```

**Winner:** Mixed-effects (uses more data efficiently)

### 3. Interpretation

**Paired Analysis:**
- Simple: "Metabolite increased by X units from pre to active"
- Clinically intuitive
- Easy to explain to non-statisticians

**Mixed-Effects:**
- More complex: "After accounting for patient-specific baselines and time effects, 2R+ status is associated with X unit increase"
- Population-level effect (average across all patients)
- Requires more statistical sophistication to interpret

**Winner:** Paired (simpler), but both have value

### 4. Assumptions

**Paired Analysis:**
- Assumes differences are normally distributed (or uses non-parametric alternative)
- No assumption about within-patient correlation (inherently paired)
- Robust to violations

**Mixed-Effects:**
- Assumes linear relationships
- Assumes residuals are normally distributed
- Assumes random effects are normally distributed
- More assumptions = more ways to be wrong, but more flexible

**Winner:** Paired (fewer assumptions)

### 5. Hypothesis Testing

**Paired Analysis (06):**
```r
# Can test:
- Is pre different from active?
- Is active different from post?
```

**Mixed-Effects (07):**
```r
# Can test:
- Is there an overall ACR effect? (controlling for time)
- Is there a time effect? (controlling for ACR)
- Does ACR affect metabolite trajectory over time? (interaction)
- Are patient-specific slopes different? (random slopes)
```

**Winner:** Mixed-effects (more flexible)

---

## Practical Examples

### Example 1: Balanced Design with 2 Timepoints

**Data:**
- 15 patients
- Each has exactly 2 samples: pre-2R and active-2R

**Recommendation:** **Use 06 (paired analysis)**

**Why:**
- Perfect for paired t-test / Wilcoxon
- No advantage to mixed-effects (same data used)
- Simpler interpretation
- More robust with small N

**Optional:** Also run 07 as confirmation (should give very similar results)

---

### Example 2: Unbalanced Design with Variable Timepoints

**Data:**
- 20 patients
- 5 patients have 2 samples (pre, active)
- 3 patients have 3 samples (pre, active, post)
- 8 patients have 4+ samples (multiple timepoints)
- 4 patients have only 1 sample

**Recommendation:** **Use 07 (mixed-effects) as PRIMARY**

**Why:**
- 06 can only use the 5+3 = 8 patients with specific paired timepoints
- 07 uses ALL 20 patients and ALL samples
- Much better power
- Can test time effects and interactions

**Optional:** Also run 06 on the 8 paired patients as sensitivity analysis

---

### Example 3: 3 Timepoints (Pre, Active, Post)

**Data:**
- 12 patients
- All have exactly 3 timepoints: pre-2R, active-2R, post-2R

**Recommendation:** **Use BOTH 06b AND 07**

**Why:**
- 06b: Gives Friedman test + reversibility patterns (clinically meaningful)
- 07: Tests whether trajectories differ (ACR×time interaction)
- Complementary information

**Interpretation:**
- 06b tells you: "8/12 patients showed full recovery"
- 07 tells you: "On average, 2R+ patients have steeper decline post-treatment (interaction p=0.03)"

---

### Example 4: Limited Sample Size

**Data:**
- 6 patients
- Variable timepoints (2-4 per patient)

**Recommendation:** **Use 06 (paired analysis)**

**Why:**
- Small N makes mixed-effects models unstable
- May not converge
- Even if converges, estimates unreliable
- Paired tests more robust with small samples

**Red Flag:** If mixed-effects model warns about convergence, trust paired analysis instead

---

## Statistical Considerations

### When Mixed-Effects Models Fail

**Symptoms:**
- "Model failed to converge"
- Extreme estimates (e.g., ACR effect = 10,000)
- Very large standard errors
- Warnings about singular fit

**Solutions:**
1. Fall back to paired analysis (06/06b)
2. Simplify model (remove random slopes, interactions)
3. Center/scale predictors
4. Check if you have enough data (need ≥10 patients ideally)

### When Paired Tests Are Underpowered

**Symptoms:**
- Few complete pairs (e.g., only 5 patients with both timepoints)
- Large variability in differences
- Non-significant results despite clear patterns

**Solutions:**
1. Try mixed-effects (07) to use all available data
2. Report effect sizes even if non-significant
3. Acknowledge limitations (small sample)

---

## Complementary Use (Recommended!)

### Best Practice: Use Both

**Why use both?**
1. **Validation:** If both methods agree, strong evidence
2. **Robustness:** Different assumptions → robust findings matter most
3. **Completeness:** Paired gives simple interpretation, mixed-effects gives sophisticated modeling
4. **Publication:** Reviewers appreciate seeing multiple approaches

**Workflow:**
```r
# Run paired analysis (simple, interpretable)
source("Scripts/06_Paired_Pre2R_vs_Active2R")

# Run mixed-effects (sophisticated, uses all data)
source("Scripts/07")
# Script 07 automatically compares results with 06

# Report:
- Metabolites significant in BOTH → Main findings (robust)
- Metabolites significant in 07 only → Hypothesis-generating (may be weak effect)
- Metabolites significant in 06 only → Strong paired effect, weaker population effect
```

### Interpretation Table

| **06 Result** | **07 Result** | **Interpretation** |
|---------------|---------------|--------------------|
| Sig | Sig | **Strong evidence** - report as main finding |
| Sig | Non-sig | Paired effect exists but weaker at population level |
| Non-sig | Sig | Population effect but high paired variability |
| Non-sig | Non-sig | No evidence of effect |

---

## Common Questions

### Q: Which is more "correct"?

**A:** Both are correct for their respective questions:
- Paired: "Do matched samples differ?"
- Mixed-effects: "Is there a population-level effect accounting for patient heterogeneity?"

Use both for comprehensive understanding.

### Q: What if they disagree?

**A:** Explore why:
1. Check sample sizes (07 using more data?)
2. Check assumptions (are residuals normal? linear relationships?)
3. Look at individual patient plots (outliers? non-linear patterns?)
4. Consider biological interpretation (is one more clinically relevant?)

Generally trust concordant findings; disagreement suggests need for caution.

### Q: Can I skip paired analysis and just use mixed-effects?

**A:** Not recommended:
- Paired tests are simpler and more interpretable
- Better for small samples
- Easier to communicate to non-statisticians
- Important for validation

Always run both if possible.

### Q: What about more complex mixed-effects models?

**A:** Script 07 fits 3 models:
1. Random intercept (simplest)
2. Random intercept + slope (more flexible)
3. Interaction model (tests ACR×time)

More complex models (e.g., crossed random effects, non-linear) possible but:
- Require larger samples
- More prone to convergence issues
- Harder to interpret
- Often unnecessary for this application

Start simple, only add complexity if justified by data and hypothesis.

---

## Decision Tree

```
START: Do you have repeated measures?
  ↓
  NO → Use 03_tests (sample-level) or 05a (patient-level)
  YES → Continue
  ↓
How many timepoints per patient?
  ↓
  Exactly 2 for all patients
    → Use 06_Paired_Pre2R_vs_Active2R (paired analysis)
    → Optionally confirm with 07 (should agree)
  ↓
  Exactly 3 for all patients
    → Use 06b_Paired_Pre_Active_Post (reversibility patterns)
    → Use 07 (trajectory modeling)
    → Report BOTH (complementary)
  ↓
  Variable (some 2, some 3, some 4+)
    → PRIMARY: Use 07 (mixed-effects, uses all data)
    → SECONDARY: Use 06 on subset with pairs (sensitivity)
  ↓
What's your sample size?
  ↓
  Small (n < 10 patients)
    → Use 06/06b only (more robust)
    → 07 may not converge
  ↓
  Moderate to large (n ≥ 10)
    → Use BOTH 06/06b AND 07
    → Compare results
    → Report concordant findings as main results
```

---

## Summary Recommendations

### ✅ Use Paired Analysis (06/06b) When:
- You have exactly 2 or 3 timepoints per patient
- Sample size is small (n < 10)
- You want simple, interpretable results
- You're interested in specific timepoint comparisons
- You want individual recovery pattern classification

### ✅ Use Mixed-Effects Models (07) When:
- You have variable numbers of timepoints per patient
- You want to use ALL available data efficiently
- Sample size is adequate (n ≥ 10 patients)
- You want to test time effects or ACR×time interactions
- You want population-level effect estimates

### ✅✅ Use BOTH When:
- You have adequate sample size
- You want robust, validated findings
- You're preparing for publication
- You want complementary perspectives

---

## Final Thoughts

**Paired analysis and mixed-effects models are not competitors - they're teammates.**

- **Paired tests** answer: "What happens within individuals?"
- **Mixed-effects** answer: "What happens at the population level, accounting for individual variation?"

Both questions are scientifically valid and clinically important. 

**Best practice:** Run both, compare results, trust concordant findings, explore discordant ones.

Your analysis pipeline is now complete with both approaches available! 🎉

---

**See also:**
- `ANALYSIS_GUIDE.md` for detailed script documentation
- `QUICK_REFERENCE.md` for quick lookup
- Script headers (06, 06b, 07) for implementation details
