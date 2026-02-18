C Engineering constants per material
      IF (MATNAME == 'EW') THEN
      KW_E1 = (0.269D0)
      KW_E2 = (0.090D0)
      KW_E3 = (9.855D0)
      KW_nu23 = (0.006D0)
      KW_nu13 = (0.018D0)
      KW_nu12 = (1.434D0)
      KW_G23 = (0.556D0)
      KW_G13 = (0.799D0)
      KW_G12 = (0.037D0)
      ELSE IF (MATNAME == 'TW') THEN
      KW_E1 = (0.614D0)
      KW_E2 = (0.416D0)
      KW_E3 = (16.493D0)
      KW_nu23 = (0.017D0)
      KW_nu13 = (0.027D0)
      KW_nu12 = (0.788D0)
      KW_G23 = (1.053D0)
      KW_G13 = (1.163D0)
      KW_G12 = (0.126D0)
      ELSE IF (MATNAME == 'LW') THEN
      KW_E1 = (1.386D0)
      KW_E2 = (2.394D0)
      KW_E3 = (29.538D0)
      KW_nu23 = (0.061D0)
      KW_nu13 = (0.028D0)
      KW_nu12 = (0.278D0)
      KW_G23 = (2.488D0)
      KW_G13 = (1.586D0)
      KW_G12 = (0.479D0)
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
