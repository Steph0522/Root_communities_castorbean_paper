make_face_names <- function(mat, rc_fun, rc_names_b = NA, 
                            rc_names_i = NA) {
  f_names <- rc_fun(mat)
  ids_b <- rc_names_b %>% match(rc_fun(mat))
  ids_i <- rc_names_i %>% match(rc_fun(mat))
  ids_bi <- rc_names_i %>% match(rc_fun(mat))
  
  ids_b %>%
    walk(
      function(i)
        f_names[i] <<-
        bquote(bold(.(rc_fun(mat)[i]))) %>%
        as.expression()
    )
  ids_i %>%
    walk(
      function(i)
        f_names[i] <<-
        bquote(italic(.(rc_fun(mat)[i]))) %>%
        as.expression()
    )
  
  f_names
}
