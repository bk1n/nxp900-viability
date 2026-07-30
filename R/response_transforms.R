ic50 <- function(intensity, negative_control) {
    (intensity / negative_control) * 100
}

gi50 <- function(intensity, negative_control, zero_control) {
    ((intensity - zero_control) / (negative_control - zero_control)) * 100
}

gr <- function(intensity, negative_control, zero_control) {
    (2^(log2(intensity / zero_control) / log2(negative_control / zero_control)) - 1) * 100
}
