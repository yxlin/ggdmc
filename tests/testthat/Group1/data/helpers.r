simple_get_accuracy <- function(df) {
    C <- ifelse(
        (df$S == "s1" & df$R == "r1") | (df$S == "s2" & df$R == "r2"),
        "O",
        "X"
    )
    prop <- table(C) / nrow(df)
    print(prop)
    return(prop)
}

complex_get_accuracy <- function(df) {
    C <- ifelse(
        (df$S == "紅" & df$R == "反應東") |
            (df$S == "黃" & df$R == "反應南") |
            (df$S == "藍" & df$R == "反應西") |
            (df$S == "綠" & df$R == "反應北"),
        "O",
        "X"
    )
    prop <- table(C) / nrow(df)
    print(prop)
    return(prop)
}

zx_get_accuracy <- function(df) {
    C <- ifelse(
        (df$S == "left" & df$R == "z_key") |
            (df$S == "right" & df$R == "x_key"),
        "O",
        "X"
    )
    prop <- table(C) / nrow(df)
    print(prop)
    return(prop)
}
