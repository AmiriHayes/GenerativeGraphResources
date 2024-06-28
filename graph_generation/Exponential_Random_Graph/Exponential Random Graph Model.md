---
affiliation: RAND Corporation
author:
- name: Raffaele Vardavas
bibliography: ../../References/ERGM.bib
date: "`r format(Sys.time(), '%B %d, %Y')`"
output: "rmdformats::readthedown"
title: Exponential family Random Graph Models (ERGMs) Approach
---

# Introduction and Aim

\`\`\`{r setup, include=FALSE, echo=FALSE} remove(list = ls())

knitr::opts_chunk\$set(cache=T, comment=NA, fig.align='center')

# Needed packages to perform the analysis

load \<- function(pkg){ new.pkg \<- pkg\[!(pkg %in%
installed.packages()\[, "Package"\])\] if (length(new.pkg))
install.packages(new.pkg, dependencies = TRUE) sapply(pkg, require,
character.only = TRUE) }

packages \<-
c("knitr","kableExtra","rmdformats","rmarkdown","dplyr","here",
"DBI","igraph","rgdal","zipcode","cowplot","DT",
"htmltools","visNetwork","reticulate",
"statnet","intergraph","ergm","ergm.ego") load(packages)

library.dir \<- "R Notebooks/library_network_analysis/" library \<-
here(library.dir,"library.R") source(library)


    This notebook explores Exponential family Random Graph Models (ERGMs) approach and how to use them to large-scale network data and merge in ego-centric data can from a different dataset. We focus on the how to implement statistical analysis of network data with ERGMs using the `statnet` software tools. This notebook copies various sections in its entirely from the existing workshops notebooks available at https://github.com/statnet/Workshops/wiki. All the work here is very much based off from those notebooks. As such this Notebook is for internal use only. For more detail please see the [references](ergm_tutorial.html#references) at the end of this tutorial.  


    # Brief Introduction to ERGMs

    Exponential-family random graph models (ERGMs) are a general class of models based in exponential-family theory for specifying the probability distribution for a set of random graphs or networks. Within this framework, one can---among other tasks:

    * Define a model for a network that includes covariates representing features like homophily, mutuality, triad effects, and a wide range of other structural features of interest;

    * Obtain maximum-likehood estimates for the parameters of the specified model for a given data set; 
    * Test individual coefficients, assess models for convergence and goodness-of-fit and perform various types of model comparison; and 

    * Simulate new networks from the underlying probability 
    distribution implied by the fitted model.

    ### The general form for an ERGM 

    ERGMs are a class of models, like linear regression or GLMs.  The general form of the model specifies the probability of the entire network (the left hand side), as a function of terms that represent network features we hypothesize may occur more or less likely than expected by chance (the right hand side).  The general form of the model can be written as:

    $$
    P(Y=y)=\frac{\exp(\theta'g(y))}{k(\theta)}
    $$

    where 

    * $Y$ is the random variable for the state of the network (with realization $y$), 

    * $g(y)$ is a vector of model statistics for network $y$, 

    * $\theta$ is the vector of coefficients for those statistics, and 
    * $k(\theta)$ represents the quantity in the numerator summed over all possible networks (typically constrained to be all networks with the same node set as $y$).

    If you're not familiar with the compact notation here, note that the numerator represents a formula that is linear in the log form:

    $$
    \log({\exp(\theta'g(y))}) = \theta_1g_1(y) + \theta_2g_2(y)+ ... + \theta_pg_p(y)
    $$
    where $p$ is the number of terms in the model.  From this one can more easily observe the analogy to a traditional statistical model:  the coefficients $\theta$ represent the size and direction of the effects of the covariates $g(y)$ on the overall probability of the network.

    ### The model statistics $g(y)$:  ERGM terms

    The statistics $g(y)$ can be thought of as the "covariates" in the model.  In the network modeling context, these represent network features like density, homophily, triads, etc.  In one sense, they are like covariates you might use in other statistical models.  But they are different in one important respect:  these $g(y)$ statistics are functions of the network itself -- each is defined by the frequency of a specific configuration of dyads observed in the network -- so they are not measured by a question you include in a survey (e.g., the income of a node), but instead need to be computed on the specific network you have, after you have collected the data.

    As a result, every term in an ERGM must have an associated algorithm for computing its value for your network.  The `ergm` package in `statnet` includes about 150 term-computing algorithms.  We will explore some of these terms in this 
    tutorial, and links to more information are provided in 
    [section 3.](ergm_tutorial.html#model-terms-available-for-ergm-estimation-and-simulation).

    You can get the list of all available terms, and the syntax for using them, by typing:

    ```{r, eval=FALSE}
    ?'ergm-terms'

You can also search for terms with keywords:

`{r eval=FALSE, include=FALSE} search.ergmTerms(keyword='homophily')`

For more information, see the vignette on `ergm-terms`:
`{r eval=FALSE, include=FALSE} vignette('ergm-term-crossRef')`

One key distinction in model terms is worth keeping in mind: terms are
either *dyad independent* or *dyad dependent*. Dyad *independent* terms
(like nodal homophily terms) imply no dependence between dyads---the
presence or absence of a tie may depend on nodal attributes, but not on
the state of other ties. Dyad *dependent* terms (like degree terms, or
triad terms), by contrast, imply dependence between dyads. Such terms
have very different effects, and much of what is different about network
models comes from these terms. They introduce complex cascading effects
that can often lead to counter-intuitive and highly non-linear outcomes.
In addition, a model with dyad dependent terms requires a different
estimation algorithm, so when we use them below you will see some
different components in the output.

### ERGM probabilities: at the tie-level

The ERGM expression for the probability of the entire graph shown above
can be re-expressed in terms of the conditional log-odds of a single tie
between two actors:

$$
\operatorname{logit}{(Y_{ij}=1|y^{c}_{ij})=\theta'\delta(y_{ij})}
$$

where

-   $Y_{ij}$ is the random variable for the state of the actor pair
    $i,j$ (with realization $y_{ij}$), and

-   $y^{c}_{ij}$ signifies the complement of $y_{ij}$, i.e. all dyads in
    the network other than $y_{ij}$.

-   $\delta(y_{ij})$ is a vector of the "change statistics" for each
    model term. The change statistic records how the $g(y)$ term changes
    if the $y_{ij}$ tie is toggled on or off. So:

$$
\delta(y_{ij}) = g(y^{+}_{ij})-g(y^{-}_{ij})
$$

where

-   $y^{+}_{ij}$ is defined as $y^{c}_{ij}$ along with $y_{ij}$ set to
    1, and
-   $y^{-}_{ij}$ is defined as $y^{c}_{ij}$ along with $y_{ij}$ set to
    0.

So $\delta(y_{ij})$ equals the value of $g(y)$ when $y_{ij}=1$ minus the
value of $g(y)$ when $y_{ij}=0$, but all other dyads are as in $y$.

This expression shows that the coefficient $\theta$ can be interpreted
as that term's contribution to the log-odds of an individual tie,
conditional on all other dyads remaining the same. The coefficient for
each term in the model is multiplied by the number of configurations
that tie will create (or remove) for that specific term. We will see
exactly how this works in the sections that follow.

# ERGMs of Ego-Centric Data

Below are important points taken from the workshop notes:

1.  One of the most powerful features of ERGMs is that they can be used
    to estimate models from from egocentrically sampled data, and the
    fitted models can then be used to simulate complete networks (of any
    size) that will have the properties of the original network that are
    observed and represented in the model.

2.  In many empirical contexts, it is not feasible to collect a network
    census or even an adaptive (link-traced) sample. Even when one of
    these may be possible in practice, egocentrically sampled data are
    typically cheaper and easier to collect.

3.  Long regarded as the poor country cousin in the network data family,
    egocentric data contain a remarkable amount of information. With the
    right statistical methods, such data can be used to explore the
    properties of the complete networks in which they are embedded. The
    basic idea here is to combine what is observed, with assumptions, to
    define a class of models that describe the distribution of networks
    that are centered on the observed properties. The variation in these
    networks quantifies some of the uncertainty introduced by the
    assumptions.

4.  The egocentric estimation/simulation framework extends to temporal
    ERGMs ("TERGMs") as well, with the minimal addition of an estimate
    of partnership duration. This makes it possible to simulate complete
    dynamic networks from a single cross-sectional egocentrically
    sampled network. For an example of what you can do with this, check
    out the network movie we developed to explore the impact of dynamic
    network structure on HIV transmission, see http://statnet.org/movies

5.  In this notebook we will use the specific package `ergm.ego` and the
    material given in the workshop that can be found online at the
    [statnet Workshops wiki](https://github.com/statnet/Workshops/wiki).

6.  It is very important to note that for the moment `ergm.ego` uses the
    minimal egocentric network study design, in which alters cannot be
    uniquely identified and alter matrices are not collected. The
    minimal design is more common, and the data are more widely
    available, largely because it is less invasive and less
    time-consuming than designs which include identifiable alter
    matrices. However, development of estimation where alter--alter
    matrices are available is being planned.

# Simulating a Social Network of a School High School

Let's use a simulated mutual friendship network based on one of the
schools from the Add-Health study. Here, we'll examine the homophily in
friendships by grade and race. Both are discrete attributes so we use
the ergm-term ***nodematch***.

``` {r}
data(faux.mesa.high) 
mesa <- faux.mesa.high
```

``` {r}
mesa
par(mfrow=c(1,1)) # Back to 1-panel plots
plot(mesa, vertex.col='Grade')
legend('bottomleft',fill=7:12,
       legend=paste('Grade',7:12),cex=0.75)
```

## A Simple ERGM

Let's proceed with an ERGM analysis that models the network based on the
total number of edges (*edges*) in the network, $\sum{y_{ij}}$. The name
of this ergm-term is `edges`, and when included in an ERGM its
coefficient controls the overall density of the network.

`{r message=FALSE, warning=FALSE} fauxmodel.00 <- ergm(mesa ~edges ) summary(fauxmodel.00)`

How should we interpret the coefficients? The easiest way is to return
to the logit form of the ERGM. The log-odds that a tie is present is $$
\small{
\begin{eqnarray*}
logit(p(y)) &=& \theta \times \delta(g(y)) 
\end{eqnarray*}
},
$$ where $\delta(g(y))$ here is the change in the covariates. So if we
assume that these are all unit values, we can find the corresponding
probability is obtained by taking inverse logit of $\theta$:

``` {r}
print(plogis(fauxmodel.00$coef) )
```

This probability corresponds to the density we observe in the network:
there are 203 ties and $205\choose{2}=20910$ dyads, so the probability
of a tie is $203/20910=0.0097$.

## A Less-Simple ERGM

Now lets add two additional covariates

1.  Factor attribute (*nodefactor*) effects on Grade and Race: This
    specifies the dependence one or more categorical attributes of the
    vertices.
2.  Differential homophily (*nodematch*): Here p network statistics are
    added to the model, where p is the number of unique values of the
    attr attribute. The kth such statistic counts the number of edges
    (i,j) for which attr(i) == attr(j) == value(k), where value(k) is
    the kth smallest unique value of the attr attribute. This is also
    called differential homophily, because each group is allowed to have
    a unique propensity for within-group ties.

`{r warning=FALSE, include=FALSE} fauxmodel.01 <- ergm(mesa ~edges +                         nodefactor('Grade') + nodematch('Grade',diff=T) +                        nodefactor('Race') + nodematch('Race',diff=T)) summary(fauxmodel.01) print(plogis(fauxmodel.01$coef) )`

Note:

1.  The coefficient for the edge has changed. The reason is that now the
    creation of edges at random is much less likely and it is controlled
    by the homophily.
2.  That two of the coefficients are estimated as -Inf (the nodematch
    coefficients for race Black and Other). Why is this?

``` {r}
table(mesa %v% 'Race') # Frequencies of race
mixingmatrix(mesa, "Race")
```

The problem is that there are very few students in the Black and Other
race categories, and these few students form no within-group ties. The
empty cells are what produce the -Inf estimates.

Note that we would have caught this earlier if we had looked at the
$g(y)$ stats at the beginning:

``` {r}
summary(mesa ~edges  + 
          nodefactor('Grade') + nodematch('Grade',diff=T) +
          nodefactor('Race') + nodematch('Race',diff=T))
```

**Moral**: It's important to check the descriptive statistics of a model
in the observed network before fitting the model.

See also the ergm-term ***nodemix*** for fitting mixing patterns other
than homophily on discrete nodal attributes.

## The Ego-Centric Network Data

Now, let's turn this into an ego-data object:
`{r warning=FALSE} mesa.ego <- as.egodata(mesa) # Generates warning because there are no vertex IDs.`

``` {r}
str(mesa.ego)
summary(mesa.ego)
head(mesa.ego$egos)
head(mesa.ego$alters)
```

Note that the ego table contains the egoID, and the nodal attributes
Race, Grade and Sex. The alter table also contains the egoID, an
equivalent set of nodal attributes. The egoID column is used by
functions that know how to work with ego-data objects to match alters to
the ego that nominated them. The alters do not have unique IDs (because
we do not know if they are unique).

We can compare the mixing matrices in the full data and the egocentric
data"

``` {r}
mixingmatrix(mesa, "Grade")
mixingmatrix(mesa.ego, "Grade")
round(100*mixingmatrix(mesa.ego, "Grade", rowprob=T),1)
```

## A Simple Ego-ERGM

Let's start with simple edges-only model to see what's the same and what
is different from a call to ergm

`{r message=FALSE, warning=FALSE} fit.edges <- ergm.ego(mesa.ego ~ edges) summary(fit.edges)`

``` {r}
names(fit.edges)
fit.edges$ppopsize
fit.edges$popsize
```

Many of the elements of the object are the same as you would get from an
ergm fit, but the last few elements are unique to ergm.ego. Here you can
see the ppopsize -- the pseudo-population size used to construct the
target statistics, and popsize -- the final scaled population size after
network size adjustment is applied. The values that were used in the fit
were the default values, since we did not specify otherwise. So,
ppopsize=205 (the sample size, or number of egos), and popsize=1, so the
scaling returns the per capita estimates from the model parameters.

As the output shows, the model fit was fit using MCMC. This, too is
different from the edges-only model using ergm. For ergm, models with
only dyad-dependent terms are fit using Newton-Raphson algorithms (the
same algorithm used for logistic regression), not MCMC. For ergm.ego,
estimation is away based on MCMC, regardless of the terms in the model.

Now let's see what the MCMC diagnostics for this model look like

``` {r}
mcmc.diagnostics(fit.edges)
```

We then do a goodness-of-fit (GOF) analysis of the overall model and by
comparing observed statistics that are not in the model, like the full
degree distribution, with simulations from the fitted model. This is the
same procedure that we use for ergm, but now with a more limited set of
observed higher-order statistics to use for assessment.

``` {r}
plot(gof(fit.edges, GOF="model"))
plot(gof(fit.edges, GOF="degree"))
```

And here, finally, we see some bad behavior, but this too is expected
from such a simple model. The GOF plot shows there are almost twice as
many isolates in the observed data than would be predicted from a simple
edges-only model. Of course we knew this from having looked at the
degree distribution plots with the Bernoulli random graph overlay.

Ok, so that's a full cycle of description, estimation, and model
assessment.

From here, let's try fitting a degree(0) term to see how that changes
the degree distribution assessment.

`{r message=FALSE, warning=FALSE} fit.deg0 <- ergm.ego(mesa.ego ~ edges + degree(0), control=control.ergm.ego(ppopsize=1000)) summary(fit.deg0)`

``` {r}
plot(gof(fit.deg0, GOF="model"))
plot(gof(fit.deg0, GOF="degree"))
```

## A Less-Simple Ego-ERGM

So, we've now fit the isolates exactly, and the fit is better, but this
now suggests there are more nodes with just one tie than would be
expected, given the mean degree, and the number of isolates.

And just to round things off, let's fit a relatively large model. Note
that we specify the omitted group for the Race nodefactor term with the
levels argument. In general, it's a good idea to omit the largest group,
so that the comparisons have a robust base. That's group "2" here,
Hispanics.

\`\`\`{r message=FALSE, warning=FALSE} set.seed(1)

fit.full \<- ergm.ego(mesa.ego \~ edges + degree(0:1) +
nodefactor("Sex") + nodefactor("Race", levels=-2) +
nodefactor("Grade") + nodematch("Sex") + nodematch("Race") +
nodematch("Grade")) summary(fit.full)

    ```{r}
    #mcmc.diagnostics(fit.full)
    plot(gof(fit.full, GOF="model"))
    plot(gof(fit.full, GOF="degree"))

## Generating Networks from the Model

It is also possible to simulate complete networks from this fit:

``` {r}
sim.full <- simulate(fit.full)

# compare values of the observed to the simulated statistics
cbind( original = summary(mesa ~ edges + degree(0:1)
                          + nodefactor("Sex") + nodefactor("Race", levels=-2) + nodefactor("Grade")
                          + nodematch("Sex") + nodematch("Race") + nodematch("Grade")),
       simulated = summary(sim.full ~ edges + degree(0:1)
                           + nodefactor("Sex") + nodefactor("Race", levels=-2) + nodefactor("Grade")
                           + nodematch("Sex") + nodematch("Race") + nodematch("Grade")),
       difference = summary(mesa ~ edges + degree(0:1)
                            + nodefactor("Sex") + nodefactor("Race", levels=-2) + nodefactor("Grade")
                            + nodematch("Sex") + nodematch("Race") + nodematch("Grade")) -
                    summary(sim.full ~ edges + degree(0:1)
                            + nodefactor("Sex") + nodefactor("Race", levels=-2) + nodefactor("Grade")
                            + nodematch("Sex") + nodematch("Race") + nodematch("Grade")))
```

``` {r}
# eyeball the plot of the simulated data
plot(sim.full, vertex.col="Grade")
legend('bottomleft',fill=7:12,legend=paste('Grade',7:12),cex=0.75)
```

## Generating Larger Networks

We can simulate a larger school by changing the ppopsize.

\`\`\`{r message=FALSE, warning=FALSE} set.seed(1)

fit.full.larger \<- ergm.ego(mesa.ego \~ edges + degree(0:1) +
nodefactor("Sex") + nodefactor("Race", levels=-2) +
nodefactor("Grade") + nodematch("Sex") + nodematch("Race") +
nodematch("Grade"), control=control.ergm.ego(ppopsize=1000))
sim.full.larger \<- simulate(fit.full.larger) \# eyeball the plot of the
simulated data plot(sim.full.larger , vertex.col="Grade")
legend('bottomleft',fill=7:12,legend=paste('Grade',7:12),cex=0.75)



    # Attempt to Network Data Fuse

    Here we begin to experiment with Ego-ERGM. This is our material, departing from the steps followed in the publicly available workshop notebooks. Here we ask if we can produce two different Ego-ERGMs on different representations of the same data and fuse the model? Let's start from our simulated school with 1000 students, and construct a new node attribute describing behavior. Let's call it the Talkative attribute. How can we generate new data that is not colinear with Grade, Sex and Race? Well perhaps by mixing in network structure. Here are the steps we follow:

    1. Generate a dummy numerical attribute that depends on Grade, Sex and Race - and hence is co-linear.
    2. Calculate it's mean over the student population.
    3. For each student calculate the average value over all his/her alters. This replaces the dummy numerical attribute.
    4. For isolates, replace their value with the mean found in step 2.
    5. Finally, compute the quantiles of the new dummy attribute and use them to provide a categorical attribute with 4 levels labeled vL, L, M and H.

    This generates our Talkative attribute. Note, a students Talkative attribute depends on his/her degree and not just on Grade, Sex and Race. Hence it is no longer co-linear. 

    ```{r message=FALSE, warning=FALSE}
    g.mesa <-asIgraph(sim.full.larger) ## this converts the mesa data to igraph format
    adj.list<-adjacent_vertices(g.mesa,V(g.mesa))

    tlk<-as.numeric(as.factor(V(g.mesa)$Race))+as.numeric(as.factor(V(g.mesa)$Sex))+as.numeric(V(g.mesa)$Grade)/12
    tlk<- tlk-min(tlk)  
    tlk <- tlk/max(tlk)
    V(g.mesa)$Talkative <- tlk

    tlk <-sapply(adj.list, FUN=function(x){
      return(mean(V(g.mesa)$Talkative[x]))
    })


    default.grade <- mean(V(g.mesa)$Grade)
    default.mean.value <- mean(tlk[!is.nan(tlk)]) 
    tlk[is.nan(tlk)] <- default.mean.value  

    cor.val <- cor(tlk, V(g.mesa)$Grade)

    default.tlk <- cut(default.mean.value,breaks=quantile(tlk, probs=c(0:4)/4),include.lowest =T)
    levels(default.tlk) <- c("vL","L","M","H")

    tlk <- cut(tlk,breaks=quantile(tlk, probs=c(0:4)/4),include.lowest =T)
    levels(tlk) <- c("vL","L","M","H")

    V(g.mesa)$Talkative <- as.character(tlk)
    mesa.new.big <- asNetwork(g.mesa)
    mesa.ego.new.big <- as.egodata(mesa.new.big)

Now let's assume that we have two different school networks. The first
is our larger school with 1000 students. But the data for the Talkative
attribute does not exist. Hence, each student in this data-set has an NA
attribute for Talkative. The second is our smaller school with 205
nodes. Here however although the Talkative attribute is available, we do
not have Grade information and these are set to NA.

\`\`\`{r message=FALSE, warning=FALSE} g.mesa \<-asIgraph(mesa) \## this
converts the mesa data to igraph format
adj.list\<-adjacent_vertices(g.mesa,V(g.mesa))

tlk\<-as.numeric(as.factor(V(g.mesa)$Race))+as.numeric(as.factor(V(g.mesa)$Sex))+as.numeric(V(g.mesa)$Grade) tlk<- tlk-min(tlk) tlk <- tlk/max(tlk) V(g.mesa)$Talkative
\<- tlk

tlk \<-sapply(adj.list, FUN=function(x){
return(mean(V(g.mesa)\$Talkative\[x\])) })

default.grade \<- mean(V(g.mesa)\$Grade) default.mean.value \<-
mean(tlk\[!is.nan(tlk)\]) tlk\[is.nan(tlk)\] \<- default.mean.value

cor.val \<- cor(tlk, V(g.mesa)\$Grade)

default.tlk \<- cut(default.mean.value,breaks=quantile(tlk,
probs=c(0:4)/4),include.lowest =T) levels(default.tlk) \<-
c("vL","L","M","H")

tlk \<- cut(tlk,breaks=quantile(tlk, probs=c(0:4)/4),include.lowest =T)
levels(tlk) \<- c("vL","L","M","H")

V(g.mesa)\$Talkative \<- as.character(tlk) mesa.new.small \<-
asNetwork(g.mesa) mesa.ego.new.small \<- as.egodata(mesa.new.small)


    Unfortunately, due to degeneracy errors that the ergm.ego would generate using a model that includes missing data with NAs we need an alternative approach for labeling the missing data. A first approach is to set the missing attribute to a single value representing the neutral value (possibly a mean or median value). This however also leads to degeneracy. Hence we resort to randomly scrambling the attribute. So in our smaller school randomly assign a Grade. 
    ```{r}
    mesa.ego.target.small<- mesa.ego.new.small
    mesa.ego.new.small$egos  <- mesa.ego.new.small$egos %>% 
      mutate(Grade=sample(unique(Grade),n(),replace = T) )
    mesa.ego.new.small$alters <- mesa.ego.new.small$alters %>% 
      mutate(Grade=sample(unique(Grade),n(),replace = T))

We scramble the Talkative attribute in the larger network in the same
way. However, before scrabbling it, we make a copy of it as this new
larger non scrambled network represents our target network that has all
four node attributes

``` {r}
mesa.ego.target.big<- mesa.ego.new.big
mesa.ego.new.big$egos <- mesa.ego.new.big$egos %>% 
  mutate(Talkative=sample(unique(Talkative),n(),replace = T) )
mesa.ego.new.big$alters <- mesa.ego.new.big$alters %>% 
  mutate(Talkative=sample(unique(Talkative),n(),replace = T) )
```

Let's create our ergm.ego models for our target networks, the small
school network with scrambled Grades, and the large network with
scrambled Talkative attribute.\
\`\`\`{r message=FALSE, warning=FALSE} fit.target \<-
ergm.ego(mesa.ego.target.big \~ edges + degree(0:1) +
nodefactor("Sex") + nodefactor("Race", levels=-2) +
nodefactor("Grade") + nodefactor("Talkative") + nodematch("Sex") +
nodematch("Race") + nodematch("Grade") + nodematch("Talkative"))

fit.new.small \<- ergm.ego(mesa.ego.new.small \~ edges + degree(0:1) +
nodefactor("Sex") + nodefactor("Race", levels=-2) +
nodefactor("Grade") + nodefactor("Talkative") + nodematch("Sex") +
nodematch("Race") + nodematch("Grade") + nodematch("Talkative"))

fit.new.big \<- ergm.ego(mesa.ego.new.big \~ edges + degree(0:1) +
nodefactor("Sex") + nodefactor("Race", levels=-2) +
nodefactor("Grade") + nodefactor("Talkative") + nodematch("Sex") +
nodematch("Race") + nodematch("Grade") + nodematch("Talkative"))

    Let's simulate our two networks and compare them
    ```{r}
    sim.new.small <- simulate(fit.new.small)
    sim.new.big <- simulate(fit.new.big)

``` {r}
# compare values of the observed to the simulated statistics
x<-cbind( target.small= summary(mesa.ego.target.small~ edges + degree(0:1)
                          + nodefactor("Sex") + nodefactor("Race", levels=-2) + nodefactor("Grade")+ nodefactor("Talkative")
                          + nodematch("Sex") + nodematch("Race") + nodematch("Grade")+ nodematch("Talkative")),
       simulated.new.small = summary(sim.new.small ~ edges + degree(0:1)
                                + nodefactor("Sex") + nodefactor("Race", levels=-2) + nodefactor("Grade")+ nodefactor("Talkative")
                                + nodematch("Sex") + nodematch("Race") + nodematch("Grade")+ nodematch("Talkative")),
       target.big= summary(mesa.ego.target.big~ edges + degree(0:1)
                          + nodefactor("Sex") + nodefactor("Race", levels=-2) + nodefactor("Grade")+ nodefactor("Talkative")
                          + nodematch("Sex") + nodematch("Race") + nodematch("Grade")+ nodematch("Talkative")),
       simulated.new.big = summary(sim.new.big ~ edges + degree(0:1)
                                + nodefactor("Sex") + nodefactor("Race", levels=-2) + nodefactor("Grade")+ nodefactor("Talkative")
                                + nodematch("Sex") + nodematch("Race") + nodematch("Grade")+ nodematch("Talkative")))


kable(x)
kable(round(x/colSums(x),2))
```

Both the small and the big networks are decently modeling their
respective targets.

## Fusing the two models

Now let's attempt building a new model based on the ermg.ego models for
the small and the large network and fuse them. Our simple approach is to
get the Talkative regression coefficients from the ergm.ego model of the
large school (i.e., the one with the missing Talkative attribute) and
the replace it with the Talkative regression coefficients from the
ergm.ego model of the smaller school (i.e., with missing Grade
attribute). Then, we will generate a school of size 1000 students and
see how it does. Our approach is loosely similar to the approach taken
by @hu_combining_2015 for fusing and comparing regression models.

``` {r}
fit.Combined <- fit.new.big

tmp<- grepl("Talkative",names(fit.Combined$coef))# | grepl("Grade",names(fit.Combined$coef)) 
fit.Combined$coef[tmp] <- fit.new.small$coef[tmp]#(fit.new.big$coef+fit.new.small$coef)[tmp]/2#(5+cor.val)
sim.Combined <- simulate(fit.Combined)
```

``` {r}
# compare values of the observed to the simulated statistics
x<-cbind(
       target.big= summary(mesa.ego.target.big~ edges + degree(0:1)
                          + nodefactor("Sex") + nodefactor("Race", levels=-2) + nodefactor("Grade")+ nodefactor("Talkative")
                          + nodematch("Sex") + nodematch("Race") + nodematch("Grade")+ nodematch("Talkative")),
       simulated.Combined = summary(sim.Combined ~ edges + degree(0:1)
                                + nodefactor("Sex") + nodefactor("Race", levels=-2) + nodefactor("Grade")+ nodefactor("Talkative")
                                + nodematch("Sex") + nodematch("Race") + nodematch("Grade")+ nodematch("Talkative")))


kable(x)
kable(round(x/colSums(x),2))
```

`{r message=FALSE, warning=FALSE} # eyeball the plot of the simulated data plot(sim.Combined , vertex.col="Talkative")`
CAN WE GENERATE AN ERGM MODEL WITH A CONSTRAINT THAT WE FIX THE VALUES
OF SOME KNOWN COVARIATES?
