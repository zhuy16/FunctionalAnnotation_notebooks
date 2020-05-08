#' Chapter 13 dplyr verbs for clusterProfiler

library(DOSE)
data(geneList)
de = names(geneList)[1:100]
x = enrichDO(de)

#'
#' 13.1 filter
#'
library(clusterProfiler.dplyr)

filter(x, p.adjust < .05, qvalue < 0.2)

#' 13.2 arrange
#'
#'
mutate(x, geneRatio = parse_ratio(GeneRatio)) %>%
  arrange(desc(geneRatio))

#' 13.3 select
#'
select(x, -geneID) %>% head

#' 13.4 mutate
#'
y <- mutate(x, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
y

library(ggplot2)
library(forcats)
library(enrichplot)

ggplot(y, showCategory = 20, 
       aes(richFactor, fct_reorder(Description, richFactor))) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  scale_size_continuous(range=c(2, 10)) +
  theme_minimal() + 
  xlab("rich factor") +
  ylab(NULL) + 
  ggtitle("Enriched Disease Ontology")

mutate(x, FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))

#' 13.5 slice

library(ReactomePA)
data(geneList)
x <- gsePathway(geneList)

library(clusterProfiler.dplyr)
y <- arrange(x, abs(NES)) %>% 
  group_by(sign(NES)) %>% 
  slice(1:5)

library(forcats)
library(ggplot2)
library(ggstance)
library(enrichplot)

ggplot(y, aes(NES, fct_reorder(Description, NES), fill=qvalues), showCategory=10) + 
  geom_barh(stat='identity') + 
  scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
  theme_minimal() + ylab(NULL)


#' 13.6 summarise
#'
pi=seq(0, 1, length.out=11)

mutate(x, pp = cut(pvalue, pi)) %>%
  group_by(pp) %>% 
  summarise(cnt = n()) %>% 
  ggplot(aes(pp, cnt)) + geom_col() + 
  theme_minimal() +
  xlab("p value intervals") +
  ylab("Frequency") + 
  ggtitle("p value distribution")
