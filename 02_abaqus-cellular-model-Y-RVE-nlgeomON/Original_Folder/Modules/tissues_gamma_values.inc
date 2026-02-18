C Combined gamma values from all tissues
        IF (MATNAME == 'EW') THEN
        KW_Gamma_i(1) = 19.40432D0
        KW_Gamma_i(2) = 3.42063D0
        KW_Gamma_i(3) = 24.31547D0
        KW_Gamma_i(4) = 2.57760D0
        ELSE IF (MATNAME == 'TW') THEN
        KW_Gamma_i(1) = 18.77668D0
        KW_Gamma_i(2) = 3.36306D0
        KW_Gamma_i(3) = 23.84677D0
        KW_Gamma_i(4) = 2.51269D0
        ELSE IF (MATNAME == 'LW') THEN
        KW_Gamma_i(1) = 18.61503D0
        KW_Gamma_i(2) = 3.32749D0
        KW_Gamma_i(3) = 23.64455D0
        KW_Gamma_i(4) = 2.46977D0
        END IF
