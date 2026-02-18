C Engineering constants per material
      IF (MATNAME == 'EW') THEN
      KW_E1 = (0.259D0)
      KW_E2 = (0.090D0)
      KW_E3 = (9.936D0)
      KW_nu23 = (0.005D0)
      KW_nu13 = (0.017D0)
      KW_nu12 = (1.411D0)
      KW_G23 = (0.556D0)
      KW_G13 = (0.799D0)
      KW_G12 = (0.035D0)
      ELSE IF (MATNAME == 'TW') THEN
      KW_E1 = (0.603D0)
      KW_E2 = (0.414D0)
      KW_E3 = (16.662D0)
      KW_nu23 = (0.016D0)
      KW_nu13 = (0.026D0)
      KW_nu12 = (0.791D0)
      KW_G23 = (1.053D0)
      KW_G13 = (1.164D0)
      KW_G12 = (0.125D0)
      ELSE IF (MATNAME == 'LW') THEN
      KW_E1 = (1.383D0)
      KW_E2 = (2.396D0)
      KW_E3 = (29.950D0)
      KW_nu23 = (0.059D0)
      KW_nu13 = (0.027D0)
      KW_nu12 = (0.278D0)
      KW_G23 = (2.489D0)
      KW_G13 = (1.592D0)
      KW_G12 = (0.482D0)
      ELSE IF (MATNAME == 'TEST') THEN
      KW_E1 = (1.0D0)
      KW_E2 = (1.0D0)
      KW_E3 = (1.0D0)
      KW_nu23 = (0.3D0)
      KW_nu13 = (0.3D0)
      KW_nu12 = (0.3D0)
      KW_G23 = (0.3846D0)
      KW_G13 = (0.3846D0)
      KW_G12 = (0.3846D0)
      ELSE
      PRINT *, 'Error: Material "', MATNAME, '" not supported.'
      RETURN
      END IF
