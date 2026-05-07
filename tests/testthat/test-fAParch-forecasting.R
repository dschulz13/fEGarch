window.zoo <- get("window.zoo", envir = asNamespace("zoo"))
rt <- window.zoo(SP500, end = "2002-12-31")
rt_train <- head(as.numeric(rt), -250)
.debug_mode_fEGarch <- TRUE

test_that("fAParch rolling forecast works as intended without refitting", {

    # Parametric models

    expect_equal({
      est <- afitgarch_spec() %>%
        fAParch(rt, n_test = 250)
      fcast <- predict_roll(est, parallel = FALSE)
      c(
        utils::head(zoo::coredata(fcast@cmeans), 5),
        zoo::coredata(fcast@cmeans)[101:105],
        utils::tail(zoo::coredata(fcast@cmeans), 5),
        utils::head(zoo::coredata(fcast@sigt), 5),
        zoo::coredata(fcast@sigt)[101:105],
        utils::tail(zoo::coredata(fcast@sigt), 5)
      )
    }, {
      c(-0.00140815583836231, -0.00140815583836231, -0.00140815583836231,
        -0.00140815583836231, -0.00140815583836231, -0.00140815583836231,
        -0.00140815583836231, -0.00140815583836231, -0.00140815583836231,
        -0.00140815583836231, -0.00140815583836231, -0.00140815583836231,
        -0.00140815583836231, -0.00140815583836231, -0.00140815583836231,
        0.00968678239272131, 0.00892151146524181, 0.0093829308509713,
        0.00995404645010778, 0.0101156873210945, 0.0148003844319179,
        0.0138122544847586, 0.0123133146094455, 0.0155224941926008,
        0.0170082760083993, 0.0140853253164921, 0.0134860733327627,
        0.0132997690197671, 0.0152352296553703, 0.0157835493782818)
    }, tolerance = 1e-02)

    expect_equal({
      est <- suppressWarnings(afiaparch_spec() %>%
        fAParch(rt, n_test = 250))
      fcast <- predict_roll(est, parallel = FALSE)
      c(
        utils::head(zoo::coredata(fcast@cmeans), 5),
        zoo::coredata(fcast@cmeans)[101:105],
        utils::tail(zoo::coredata(fcast@cmeans), 5),
        utils::head(zoo::coredata(fcast@sigt), 5),
        zoo::coredata(fcast@sigt)[101:105],
        utils::tail(zoo::coredata(fcast@sigt), 5)
      )
    }, {
      c(-0.0012, -0.0012, -0.0012, -0.0012, -0.0012, -0.0012, -0.0012,
        -0.0012, -0.0012, -0.0012, -0.0012, -0.0012, -0.0012, -0.0012,
        -0.0012, 0.00967, 0.00904, 0.0093, 0.00948, 0.00946, 0.01316,
        0.01225, 0.01116, 0.01545, 0.01667, 0.01304, 0.01249, 0.01212,
        0.01412, 0.01451)
    }, tolerance = 1e-02)

    expect_equal({
      est <- afigjrgarch_spec() %>%
        fAParch(rt, n_test = 250)
      fcast <- predict_roll(est, parallel = FALSE)
      c(
        utils::head(zoo::coredata(fcast@cmeans), 5),
        zoo::coredata(fcast@cmeans)[101:105],
        utils::tail(zoo::coredata(fcast@cmeans), 5),
        utils::head(zoo::coredata(fcast@sigt), 5),
        zoo::coredata(fcast@sigt)[101:105],
        utils::tail(zoo::coredata(fcast@sigt), 5)
      )
    }, {
      c(-0.00107, -0.00107, -0.00107, -0.00107, -0.00107, -0.00107,
        -0.00107, -0.00107, -0.00107, -0.00107, -0.00107, -0.00107,
        -0.00107, -0.00107, -0.00107, 0.00964, 0.00915, 0.00925, 0.00922,
        0.00912, 0.01191, 0.01114, 0.01037, 0.01574, 0.01649, 0.01221,
        0.01173, 0.01135, 0.01338, 0.01357)
    }, tolerance = 1e-02)

    expect_equal({
      est <- afingarch_spec() %>%
        fAParch(rt, n_test = 250)
      fcast <- predict_roll(est, parallel = FALSE)
      c(
        utils::head(zoo::coredata(fcast@cmeans), 5),
        zoo::coredata(fcast@cmeans)[101:105],
        utils::tail(zoo::coredata(fcast@cmeans), 5),
        utils::head(zoo::coredata(fcast@sigt), 5),
        zoo::coredata(fcast@sigt)[101:105],
        utils::tail(zoo::coredata(fcast@sigt), 5)
      )
    }, {
      c(-0.00011, -0.00011, -0.00011, -0.00011, -0.00011, -0.00011,
        -0.00011, -0.00011, -0.00011, -0.00011, -0.00011, -0.00011,
        -0.00011, -0.00011, -0.00011, 0.00996, 0.00984, 0.00966, 0.0094,
        0.0092, 0.01126, 0.01072, 0.01026, 0.01416, 0.01424, 0.01282,
        0.01215, 0.0116, 0.01247, 0.01233)
    }, tolerance = 1e-02)



    expect_equal({
      est <- fitgarch_spec() %>%
        fAParch(rt, n_test = 250)
      fcast <- predict_roll(est, parallel = FALSE)
      c(
        utils::head(zoo::coredata(fcast@cmeans), 5),
        zoo::coredata(fcast@cmeans)[101:105],
        utils::tail(zoo::coredata(fcast@cmeans), 5),
        utils::head(zoo::coredata(fcast@sigt), 5),
        zoo::coredata(fcast@sigt)[101:105],
        utils::tail(zoo::coredata(fcast@sigt), 5)
      )
    }, {
      c(-0.00126, -0.00126, -0.00126, -0.00126, -0.00126, -0.00126,
        -0.00126, -0.00126, -0.00126, -0.00126, -0.00126, -0.00126,
        -0.00126, -0.00126, -0.00126, 0.00981, 0.00907, 0.00955, 0.0101,
        0.0103, 0.01464, 0.0137, 0.01224, 0.01542, 0.01663, 0.01359, 0.01303,
        0.01285, 0.0148, 0.01521)
    }, tolerance = 1e-02)

    expect_equal({
      est <- suppressWarnings(fiaparch_spec() %>%
        fAParch(rt, n_test = 250))
      fcast <- predict_roll(est, parallel = FALSE)
      c(
        utils::head(zoo::coredata(fcast@cmeans), 5),
        zoo::coredata(fcast@cmeans)[101:105],
        utils::tail(zoo::coredata(fcast@cmeans), 5),
        utils::head(zoo::coredata(fcast@sigt), 5),
        zoo::coredata(fcast@sigt)[101:105],
        utils::tail(zoo::coredata(fcast@sigt), 5)
      )
    }, {
      c(-0.00133, -0.00133, -0.00133, -0.00133, -0.00133, -0.00133,
        -0.00133, -0.00133, -0.00133, -0.00133, -0.00133, -0.00133,
        -0.00133, -0.00133, -0.00133, 0.00986, 0.00914, 0.00957, 0.01004,
        0.01018, 0.01445, 0.0135, 0.01211, 0.01546, 0.01675, 0.01359, 0.013,
        0.01277, 0.01471, 0.01517)
    }, tolerance = 1e-02)

    expect_equal({
      est <- figjrgarch_spec() %>%
        fAParch(rt, n_test = 250)
      fcast <- predict_roll(est, parallel = FALSE)
      c(
        utils::head(zoo::coredata(fcast@cmeans), 5),
        zoo::coredata(fcast@cmeans)[101:105],
        utils::tail(zoo::coredata(fcast@cmeans), 5),
        utils::head(zoo::coredata(fcast@sigt), 5),
        zoo::coredata(fcast@sigt)[101:105],
        utils::tail(zoo::coredata(fcast@sigt), 5)
      )
    }, {
      c(-0.00128, -0.00128, -0.00128, -0.00128, -0.00128, -0.00128,
        -0.00128, -0.00128, -0.00128, -0.00128, -0.00128, -0.00128,
        -0.00128, -0.00128, -0.00128, 0.00949, 0.00905, 0.00914, 0.00902,
        0.00888, 0.01173, 0.01094, 0.01021, 0.01586, 0.01625, 0.01306,
        0.01261, 0.01221, 0.01416, 0.01416)
    }, tolerance = 1e-02)

    expect_equal({
      est <- fingarch_spec() %>%
        fAParch(rt, n_test = 250)
      fcast <- predict_roll(est, parallel = FALSE)
      c(
        utils::head(zoo::coredata(fcast@cmeans), 5),
        zoo::coredata(fcast@cmeans)[101:105],
        utils::tail(zoo::coredata(fcast@cmeans), 5),
        utils::head(zoo::coredata(fcast@sigt), 5),
        zoo::coredata(fcast@sigt)[101:105],
        utils::tail(zoo::coredata(fcast@sigt), 5)
      )
    }, {
      c(-0.00024, -0.00024, -0.00024, -0.00024, -0.00024, -0.00024,
        -0.00024, -0.00024, -0.00024, -0.00024, -0.00024, -0.00024,
        -0.00024, -0.00024, -0.00024, 0.01047, 0.0105, 0.01036, 0.00996,
        0.0096, 0.0124, 0.01158, 0.01076, 0.01281, 0.01318, 0.01378, 0.01285,
        0.012, 0.01259, 0.01263)
    }, tolerance = 1e-02)

    expect_equal({
      est <- figarch_spec() %>%
        fAParch(rt, n_test = 250)
      fcast <- predict_roll(est, parallel = FALSE)
      c(
        utils::head(zoo::coredata(fcast@cmeans), 5),
        zoo::coredata(fcast@cmeans)[101:105],
        utils::tail(zoo::coredata(fcast@cmeans), 5),
        utils::head(zoo::coredata(fcast@sigt), 5),
        zoo::coredata(fcast@sigt)[101:105],
        utils::tail(zoo::coredata(fcast@sigt), 5)
      )
    }, {
      c(-0.00015, -0.00015, -0.00015, -0.00015, -0.00015, -0.00015, -0.00015,
        -0.00015, -0.00015, -0.00015, -0.00015, -0.00015, -0.00015, -0.00015,
        -0.00015, 0.01024, 0.01018, 0.00998, 0.00963, 0.00933, 0.01184,
        0.01114, 0.0105, 0.01338, 0.0139, 0.01364, 0.01282, 0.01216, 0.01291,
        0.01297)
    }, tolerance = 1e-02)

    expect_equal({
      est <- garch_spec() %>%
        fAParch(rt, n_test = 250)
      fcast <- predict_roll(est, parallel = FALSE)
      c(
        utils::head(zoo::coredata(fcast@cmeans), 5),
        zoo::coredata(fcast@cmeans)[101:105],
        utils::tail(zoo::coredata(fcast@cmeans), 5),
        utils::head(zoo::coredata(fcast@sigt), 5),
        zoo::coredata(fcast@sigt)[101:105],
        utils::tail(zoo::coredata(fcast@sigt), 5)
      )
    }, {
      c(-0.00023, -0.00023, -0.00023, -0.00023, -0.00023, -0.00023,
        -0.00023, -0.00023, -0.00023, -0.00023, -0.00023, -0.00023,
        -0.00023, -0.00023, -0.00023, 0.01037, 0.01026, 0.01016,
        0.0099, 0.00975, 0.01198, 0.01141, 0.01093, 0.01362, 0.01279, 0.01252,
        0.01198, 0.01143, 0.01224, 0.01172)
    }, tolerance = 1e-02)


    ########


    # Semiparametric models

    expect_equal({
      est <- afitgarch_spec() %>%
        fAParch(rt, n_test = 250, use_nonpar = TRUE)
      fcast <- predict_roll(est, parallel = FALSE)
      c(
        utils::head(zoo::coredata(fcast@cmeans), 5),
        zoo::coredata(fcast@cmeans)[101:105],
        utils::tail(zoo::coredata(fcast@cmeans), 5),
        utils::head(zoo::coredata(fcast@sigt), 5),
        zoo::coredata(fcast@sigt)[101:105],
        utils::tail(zoo::coredata(fcast@sigt), 5)
      )
    }, {
      c(rep(mean(rt_train, 15)), 0.00782, 0.00723, 0.0077, 0.00834, 0.00858, 0.01295, 0.01226,
        0.01105, 0.01345, 0.01457, 0.01275, 0.01235, 0.01235, 0.01388, 0.01422)
    }, tolerance = 1e-02)

    expect_equal({
      est <- suppressWarnings(afiaparch_spec() %>%
        fAParch(rt, n_test = 250, use_nonpar = TRUE))
      fcast <- predict_roll(est, parallel = FALSE)
      c(
        utils::head(zoo::coredata(fcast@cmeans), 5),
        zoo::coredata(fcast@cmeans)[101:105],
        utils::tail(zoo::coredata(fcast@cmeans), 5),
        utils::head(zoo::coredata(fcast@sigt), 5),
        zoo::coredata(fcast@sigt)[101:105],
        utils::tail(zoo::coredata(fcast@sigt), 5)
      )
    }, {
      c(rep(mean(rt_train, 15)), 0.0079, 0.00733, 0.00774, 0.00814, 0.00825,
        0.01236, 0.01161, 0.01054, 0.01367, 0.01474, 0.01251, 0.0121, 0.01192,
        0.01354, 0.01385)

    }, tolerance = 1e-02)

    expect_equal({
      est <- afigjrgarch_spec() %>%
        fAParch(rt, n_test = 250, use_nonpar = TRUE)
      fcast <- predict_roll(est, parallel = FALSE)
      c(
        utils::head(zoo::coredata(fcast@cmeans), 5),
        zoo::coredata(fcast@cmeans)[101:105],
        utils::tail(zoo::coredata(fcast@cmeans), 5),
        utils::head(zoo::coredata(fcast@sigt), 5),
        zoo::coredata(fcast@sigt)[101:105],
        utils::tail(zoo::coredata(fcast@sigt), 5)
      )
    }, {
      c(-0.00044, -0.00044, -0.00044, -0.00044, -0.00044, -0.00044,
        -0.00044, -0.00044, -0.00044, -0.00044, -0.00044, -0.00044,
        -0.00044, -0.00044, -0.00044, 0.00798, 0.00751, 0.00774, 0.00777,
        0.00771, 0.01109, 0.01037, 0.00965, 0.0146, 0.01509, 0.01193, 0.01156,
        0.01125, 0.01306, 0.01311)
    }, tolerance = 1e-02)

    expect_equal({
      est <- afingarch_spec() %>%
        fAParch(rt, n_test = 250, use_nonpar = TRUE)
      fcast <- predict_roll(est, parallel = FALSE)
      c(
        utils::head(zoo::coredata(fcast@cmeans), 5),
        zoo::coredata(fcast@cmeans)[101:105],
        utils::tail(zoo::coredata(fcast@cmeans), 5),
        utils::head(zoo::coredata(fcast@sigt), 5),
        zoo::coredata(fcast@sigt)[101:105],
        utils::tail(zoo::coredata(fcast@sigt), 5)
      )
    }, {
      c(-0.00044, -0.00044, -0.00044, -0.00044, -0.00044, -0.00044,
        -0.00044, -0.00044, -0.00044, -0.00044, -0.00044, -0.00044, -0.00044,
        -0.00044, -0.00044, 0.00872, 0.00873, 0.00861, 0.00845, 0.0083, 0.0105,
        0.01016, 0.00986, 0.01132, 0.01175, 0.01281, 0.01241, 0.0121, 0.01235,
        0.01238)
    }, tolerance = 1e-02)



    expect_equal({
      est <- fitgarch_spec() %>%
        fAParch(rt, n_test = 250, use_nonpar = TRUE)
      fcast <- predict_roll(est, parallel = FALSE)
      c(
        utils::head(zoo::coredata(fcast@cmeans), 5),
        zoo::coredata(fcast@cmeans)[101:105],
        utils::tail(zoo::coredata(fcast@cmeans), 5),
        utils::head(zoo::coredata(fcast@sigt), 5),
        zoo::coredata(fcast@sigt)[101:105],
        utils::tail(zoo::coredata(fcast@sigt), 5)
      )
    }, {
      c(-0.00044, -0.00044, -0.00044, -0.00044, -0.00044, -0.00044, -0.00044,
        -0.00044, -0.00044, -0.00044, -0.00044, -0.00044, -0.00044, -0.00044,
        -0.00044, 0.00776, 0.00715, 0.00763, 0.00826, 0.00852, 0.0127,
        0.01201, 0.01078, 0.01319, 0.01428, 0.01215, 0.01172, 0.01171,
        0.01328, 0.01362)
    }, tolerance = 1e-02)

    expect_equal({
      est <- suppressWarnings(fiaparch_spec() %>%
        fAParch(rt, n_test = 250, use_nonpar = TRUE))
      fcast <- predict_roll(est, parallel = FALSE)
      c(
        utils::head(zoo::coredata(fcast@cmeans), 5),
        zoo::coredata(fcast@cmeans)[101:105],
        utils::tail(zoo::coredata(fcast@cmeans), 5),
        utils::head(zoo::coredata(fcast@sigt), 5),
        zoo::coredata(fcast@sigt)[101:105],
        utils::tail(zoo::coredata(fcast@sigt), 5)
      )
    }, {
      c(-0.00044, -0.00044, -0.00044, -0.00044, -0.00044, -0.00044,
        -0.00044, -0.00044, -0.00044, -0.00044, -0.00044, -0.00044,
        -0.00044, -0.00044, -0.00044, 0.00758, 0.00695, 0.0075,
        0.00843, 0.0089, 0.0131, 0.01259, 0.01122, 0.01285, 0.01379,
        0.01198, 0.01157, 0.01186, 0.01335, 0.01361)
    }, tolerance = 1e-02)

    expect_equal({
      est <- figjrgarch_spec() %>%
        fAParch(rt, n_test = 250, use_nonpar = TRUE)
      fcast <- predict_roll(est, parallel = FALSE)
      c(
        utils::head(zoo::coredata(fcast@cmeans), 5),
        zoo::coredata(fcast@cmeans)[101:105],
        utils::tail(zoo::coredata(fcast@cmeans), 5),
        utils::head(zoo::coredata(fcast@sigt), 5),
        zoo::coredata(fcast@sigt)[101:105],
        utils::tail(zoo::coredata(fcast@sigt), 5)
      )
    }, {
      c(-0.00044, -0.00044, -0.00044, -0.00044, -0.00044, -0.00044,
        -0.00044, -0.00044, -0.00044, -0.00044, -0.00044, -0.00044,
        -0.00044, -0.00044, -0.00044, 0.00801, 0.00761, 0.00781,
        0.00779, 0.00772, 0.01127, 0.0106, 0.00995, 0.01458, 0.0149,
        0.013, 0.01266, 0.01235, 0.0139, 0.01387)
    }, tolerance = 1e-02)

    expect_equal({
      est <- fingarch_spec() %>%
        fAParch(rt, n_test = 250, use_nonpar = TRUE)
      fcast <- predict_roll(est, parallel = FALSE)
      c(
        utils::head(zoo::coredata(fcast@cmeans), 5),
        zoo::coredata(fcast@cmeans)[101:105],
        utils::tail(zoo::coredata(fcast@cmeans), 5),
        utils::head(zoo::coredata(fcast@sigt), 5),
        zoo::coredata(fcast@sigt)[101:105],
        utils::tail(zoo::coredata(fcast@sigt), 5)
      )
    }, {
      c(-0.00044, -0.00044, -0.00044, -0.00044, -0.00044, -0.00044,
        -0.00044, -0.00044, -0.00044, -0.00044, -0.00044, -0.00044,
        -0.00044, -0.00044, -0.00044, 0.00895, 0.00909, 0.00897,
        0.00871, 0.0084, 0.01104, 0.01047, 0.00987, 0.0107, 0.01167,
        0.01289, 0.01212, 0.01158, 0.01171, 0.01211)
    }, tolerance = 1e-02)

    expect_equal({
      est <- figarch_spec() %>%
        fAParch(rt, n_test = 250, use_nonpar = TRUE)
      fcast <- predict_roll(est, parallel = FALSE)
      c(
        utils::head(zoo::coredata(fcast@cmeans), 5),
        zoo::coredata(fcast@cmeans)[101:105],
        utils::tail(zoo::coredata(fcast@cmeans), 5),
        utils::head(zoo::coredata(fcast@sigt), 5),
        zoo::coredata(fcast@sigt)[101:105],
        utils::tail(zoo::coredata(fcast@sigt), 5)
      )
    }, {
      c(-0.00044, -0.00044, -0.00044, -0.00044, -0.00044, -0.00044,
        -0.00044, -0.00044, -0.00044, -0.00044, -0.00044, -0.00044,
        -0.00044, -0.00044, -0.00044, 0.0088, 0.00883, 0.0087, 0.00851,
        0.00834, 0.01062, 0.01023, 0.00989, 0.01121, 0.01175, 0.01287,
        0.0124, 0.01205, 0.0123, 0.01239)
    }, tolerance = 1e-02)

    expect_equal({
      est <- garch_spec() %>%
        fAParch(rt, n_test = 250, use_nonpar = TRUE)
      fcast <- predict_roll(est, parallel = FALSE)
      c(
        utils::head(zoo::coredata(fcast@cmeans), 5),
        zoo::coredata(fcast@cmeans)[101:105],
        utils::tail(zoo::coredata(fcast@cmeans), 5),
        utils::head(zoo::coredata(fcast@sigt), 5),
        zoo::coredata(fcast@sigt)[101:105],
        utils::tail(zoo::coredata(fcast@sigt), 5)
      )
    }, {
      c(-0.00044, -0.00044, -0.00044, -0.00044, -0.00044, -0.00044,
        -0.00044, -0.00044, -0.00044, -0.00044, -0.00044, -0.00044,
        -0.00044, -0.00044, -0.00044, 0.0086, 0.00854, 0.00846, 0.00828,
        0.00817, 0.00958, 0.00917, 0.00884, 0.01092, 0.01027, 0.01002,
        0.00962, 0.00921, 0.00985, 0.00947)
    }, tolerance = 1e-02)

})
