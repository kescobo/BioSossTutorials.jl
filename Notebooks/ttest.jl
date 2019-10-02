# # T-test
#
# We often want to compare a treatment group to a control group
# to identify whether or not some response occured due to the treatment.
# An easy example is to compare the bodyweights of mice
# fed either a normal diet or a high fat diet.
#
# This example is taken from [Data Analysis for Life Scientists](https://leanpub.com/dataanalysisforthelifesciences)
#
# Control mouse population: https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleControlsPopulation.csv
# Experiment: https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleMiceWeights.csv

using CSV
using DataFrames
using Statistics
using HypothesisTests
using StatsPlots
using HTTP

controls_url = "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleControlsPopulation.csv"
hfd_url = "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleMiceWeights.csv"

controls = CSV.read(HTTP.get(controls_url).body)
hfd = CSV.read(HTTP.get(hfd_url).body)

## make a histogram of control mice
histogram(controls.Bodyweight)

#-

describe(controls.Bodyweight)

# This shows what the weights of a population of normal femal mice
# would look like, absent any special treatment.
# Now let's look at an experiemtn with 24 mice,
# where 12 mice were fed normal chow,
# and the other 12 were given a high-fat diet.

chow = hfd[hfd.Diet .== "chow", :Bodyweight]
hf = hfd[hfd.Diet .== "hf", :Bodyweight]
chowc = RGBA(0.9,0.2,0.3, 0.5)
hfc = RGBA(0.2,0.7,0.3, 0.5)

histogram(chow, color=chowc, label="chow")
histogram!(hf, color=hfc, label="high fat")

#-

boxplot([chow, hf], color=[chowc hfc], title="Bodyweights", legend=false,
         ylabel="weight (grams)", xticks=([1,2], ["chow", "high fat"]))

#-

@info "Chow"
describe(chow)
@info "High fat"
describe(hf)

# The mean of the high-fat group is clearly higher
# than the group fed normal chow.
# The question is: is this difference meaningful?
#
# Afterall, the weights substantilly overlap.
# In addition, the heaviest mouse in the highfat diet group
# is lighter than the heaviest mouse in the control population:

any(x -> x > maximum(hf), controls.Bodyweight)

# ... and about 25% of mice in the control population are heavier
# than the median high-fat diet-fed mouse.

count(x-> x > quantile(hf, .5), controls.Bodyweight) / nrow(controls)

# Of course, we usually don't have this large control population to compare to,
# but it's possible we just happened to get some particularly heavy mice
# in our high-fat group,
# or even some particularly light mice in our experimental control group.
#
# A t-test can be used to test the null hypothesis (h₀) that the two groups
# are drawn from the same population.

EqualVarianceTTest(hf, chow)
