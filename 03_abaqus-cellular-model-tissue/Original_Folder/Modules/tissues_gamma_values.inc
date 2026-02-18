C Combined gamma values from all tissues
        IF (MATNAME == 'EW') THEN
        KW_Gamma_i(1) = 19.51579D0
        KW_Gamma_i(2) = 3.51130D0
        KW_Gamma_i(3) = 25.20033D0
        KW_Gamma_i(4) = 2.69849D0
        ELSE IF (MATNAME == 'TW') THEN
        KW_Gamma_i(1) = 19.35885D0
        KW_Gamma_i(2) = 3.44311D0
        KW_Gamma_i(3) = 24.78601D0
        KW_Gamma_i(4) = 2.60936D0
        ELSE IF (MATNAME == 'LW') THEN
        KW_Gamma_i(1) = 19.24423D0
        KW_Gamma_i(2) = 3.36512D0
        KW_Gamma_i(3) = 24.15964D0
        KW_Gamma_i(4) = 2.52853D0
        END IF
