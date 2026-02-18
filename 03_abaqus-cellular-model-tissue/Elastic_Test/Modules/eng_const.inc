C Engineering constants per material
      IF (MATNAME == 'EW') THEN
      KW_E1 = (0.357D0)
      KW_E2 = (0.322D0)
      KW_E3 = (7.947D0)
      KW_nu23 = (0.027D0)
      KW_nu13 = (0.032D0)
      KW_nu12 = (0.585D0)
      KW_G23 = (0.668D0)
      KW_G13 = (0.678D0)
      KW_G12 = (0.010D0)
      ELSE IF (MATNAME == 'TW') THEN
      KW_E1 = (0.700D0)
      KW_E2 = (1.036D0)
      KW_E3 = (16.490D0)
      KW_nu23 = (0.043D0)
      KW_nu13 = (0.030D0)
      KW_nu12 = (0.334D0)
      KW_G23 = (1.305D0)
      KW_G13 = (1.270D0)
      KW_G12 = (0.071D0)
      ELSE IF (MATNAME == 'LW') THEN
      KW_E1 = (1.619D0)
      KW_E2 = (2.111D0)
      KW_E3 = (28.343D0)
      KW_nu23 = (0.053D0)
      KW_nu13 = (0.040D0)
      KW_nu12 = (0.258D0)
      KW_G23 = (2.157D0)
      KW_G13 = (2.078D0)
      KW_G12 = (0.233D0)
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
