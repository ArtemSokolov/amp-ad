## Bi-variate distribution for the Braak - CDR association
## Generated from the MSBB data using the following tidyverse query
##
## read_tsv( "msbb-wrangled.tsv.gz" ) %>% na.omit %>%
##   group_by( Braak, CDR ) %>% summarize( nSamples = n() )
##
## by Artem Sokolov

distBraakCDR <- tibble::frame_data(
                            ~Braak, ~CDR, ~nSamples,
                            0,      0.0,  14,
                            0,      0.5,   6,
                            0,      1.0,   3,
                            0,      2.0,   3,
                            0,      3.0,   4,
                            1,      0.0,  24,
                            1,      0.5,  38,
                            1,      1.0,  13,
                            1,      2.0,   3,
                            1,      3.0,   4,
                            1,      4.0,   5,
                            1,      5.0,   4,
                            2,      0.0,  25,
                            2,      0.5,  43,
                            2,      1.0,  21,
                            2,      2.0,   9,
                            2,      3.0,  14,
                            2,      4.0,   9,
                            2,      5.0,   7,
                            3,      0.0,  45,
                            3,      0.5,  23,
                            3,      1.0,  30,
                            3,      2.0,  30,
                            3,      3.0,  32,
                            3,      4.0,  11,
                            4,      0.0,  12,
                            4,      0.5,  10,
                            4,      1.0,   8,
                            4,      2.0,  17,
                            4,      3.0,  17,
                            4,      4.0,  21,
                            4,      5.0,  10,
                            5,      0.0,   3,
                            5,      0.5,   9,
                            5,      1.0,  14,
                            5,      2.0,   2,
                            5,      3.0,  35,
                            5,      5.0,  28,
                            6,      0.5,   4,
                            6,      1.0,  10,
                            6,      2.0,  67,
                            6,      3.0,  76,
                            6,      4.0,  24,
                            6,      5.0,  48 )
