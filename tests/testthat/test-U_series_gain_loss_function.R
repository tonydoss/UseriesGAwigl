result <- U_series_gain_loss_function (read.csv(system.file("extdata", "data.csv", package = "UseriesGAwigl")), nbit = 100,
                                       logT_min = 1, logT_max = 7, Tmult_min = 0.5,
                                       logk238_min = -7, logk238_max = -5, k48_min = 0.1, k48_max = 10,
                                       logf238_k238_min = -20, logf238_k238_max = 3, f48_min = 0.1, f48_max = 10)

test_that("Model returns output", {
  expect_output(str(result), "List of 7")
})
