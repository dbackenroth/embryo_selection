parental_score <- rep(runif(1000, min = 2e-04, max = 6e-04), each = 20)
var <- 3e-09
child_scores <- rnorm(length(parental_score * 20), mean = parental_score, sd = sqrt(var))

d <- data.frame(p = parental_score, c = child_scores) %>%
  mutate(fam = rep(1:1000, each = 20))

stats <- group_by(d, fam) %>%
  summarise(p = p[1], 
            v = var(c))