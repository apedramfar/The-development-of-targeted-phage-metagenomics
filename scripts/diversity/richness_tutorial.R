library("phyloseq"); packageVersion("phyloseq")

data("GlobalPatterns")
library("ggplot2"); packageVersion("ggplot2")

theme_set(theme_bw())
pal = "Set1"
scale_colour_discrete <-  function(palname=pal, ...){
  scale_colour_brewer(palette=palname, ...)
}
scale_fill_discrete <-  function(palname=pal, ...){
  scale_fill_brewer(palette=palname, ...)
}

GP <- prune_species(speciesSums(GlobalPatterns) > 0, GlobalPatterns)

plot_richness(GP)

plot_richness(GP, measures=c("Chao1", "Shannon"))

sampleData(GP)$human <- getVariable(GP, "SampleType") %in% c("Feces", "Mock", "Skin", "Tongue")

plot_richness(GP, x="human", color="SampleType", measures=c("Chao1", "Shannon"))


GPst = merge_samples(GP, "SampleType")
# repair variables that were damaged during merge (coerced to numeric)
sample_data(GPst)$SampleType <- factor(sample_names(GPst))
sample_data(GPst)$human <- as.logical(sample_data(GPst)$human)

p = plot_richness(GPst, x="human", color="SampleType", measures=c("Chao1", "Shannon"))
p + geom_point(size=5, alpha=0.7)

p$layers

p$layers <- p$layers[-1]
p + geom_point(size=5, alpha=0.7)