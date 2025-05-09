library(MASS)
library(Hotelling)

# All data verified by Andrew

# Haltica oleracea
species1 <- matrix(c(
  189, 245, 137, 163,
  192, 260, 132, 217,
  217, 276, 141, 192,
  221, 299, 142, 213,
  171, 239, 128, 158,
  192, 262, 147, 173,
  213, 278, 136, 201,
  192, 255, 128, 185,
  170, 244, 128, 192,
  201, 276, 146, 186,
  195, 242, 128, 192,
  205, 263, 147, 192,
  180, 252, 121, 167,
  192, 283, 138, 183,
  200, 294, 138, 188,
  192, 277, 150, 177,
  200, 287, 136, 173,
  181, 255, 146, 183,
  192, 287, 141, 198,
), ncol = 4, byrow = TRUE)

# Haltica carduorum
species2 <- matrix(c(
  181, 305, 184, 209,
  158, 237, 133, 188,
  184, 300, 166, 231,
  171, 273, 162, 213,
  181, 297, 163, 224,
  181, 308, 160, 223,
  177, 301, 166, 221,
  198, 308, 141, 197,
  180, 286, 146, 214,
  177, 299, 171, 192,
  176, 317, 166, 213,
  192, 312, 166, 209,
  176, 285, 141, 200,
  169, 287, 162, 214,
  164, 265, 147, 192,
  181, 308, 157, 204,
  192, 276, 154, 209,
  181, 278, 149, 235,
  175, 271, 140, 192,
  197, 303, 170, 205
), ncol=4, byrow=TRUE)

# Create data frame and group labels
data <- rbind(species1, species2)
group <- factor(c(rep("oleracea", 20), rep("carduorum", 20)))
df <- data.frame(group, data)

# Apply Hotelling's T-squared test
result <- hotelling.test(data[group == "oleracea", ], data[group == "carduorum", ])
print(result)
