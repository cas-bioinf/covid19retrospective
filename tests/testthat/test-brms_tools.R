test_that("normalize_stancode", {
  expect_equal(
    normalize_stancode("// a\nb;\n  b + c = 4; // kde\ndata"),
    normalize_stancode("// dasflkjldl\n   // adsfadsfa\n b;\n\n  \n  \t\rb + c = 4;\ndata")
  )
  expect_equal(
    normalize_stancode("data /* adfa */ {\nint a;\n /* asdddede \n asdfas \n asf */}\n"),
    normalize_stancode("data {\nint a;\n} /* aa \n adfasdf \n asdfadsf ddd */\n")
  )
  expect_equal(
    normalize_stancode("data \n {\nint a;\n\n }  \t\n"),
    normalize_stancode("data {\nint a;\n} \n")
  )
  expect_equal(
    normalize_stancode("/* \n\n */\na*/"),
    normalize_stancode("a*/")
  )
  expect_equal(
    normalize_stancode("//adsfadf \ra // asdfasdf\r\n"),
    normalize_stancode("a")
  )
  expect_equal(
    normalize_stancode("/* * \n * \n * fg / */hhh"),
    normalize_stancode("hhh")
  )
  expect_equal(
    normalize_stancode("a //b"),
    normalize_stancode("a")
  )
  expect_false(normalize_stancode("// a\ndata {\nint a;\n}\n") ==
                 normalize_stancode("// a\ndata {\nint b;\n}\n"))
  #Should not remove single whitespace
  expect_false(normalize_stancode("da ta") ==
                 normalize_stancode("data"))
  #Should handle wrong nested comments
  expect_false(normalize_stancode("/* \n\n */\na*/") ==
                 normalize_stancode("b*/"))
})
