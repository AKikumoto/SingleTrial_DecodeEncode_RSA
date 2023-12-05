# Getting probabilities in classification for each observation
# 
# =============================================================
# Load required libraries
library(caret)

# Load dataset
data("iris")

# Add noise to data
set.seed(123)
noise_factor <- 3  # Adjust the noise level as needed

iris_noisy <- iris
iris_noisy$Sepal.Length <- iris_noisy$Sepal.Length + noise_factor * rnorm(nrow(iris_noisy))
iris_noisy$Sepal.Width <- iris_noisy$Sepal.Width + noise_factor * rnorm(nrow(iris_noisy))

# Split dataset
index <- createDataPartition(iris_noisy$Species, p = 0.8, list = FALSE)
train_data <- iris_noisy[index, ]
test_data <- iris_noisy[-index, ]

# Train model and make predictions with probability estimates
model <- train(Species ~ ., data = train_data, method = "pda")
class_probabilities <- as.data.frame(predict(model, newdata = test_data, type = "prob"))

# Print or use class probabilities as needed
print(class_probabilities)

# - this class probabilities (typically way noisier) is regressed by RSAs 
