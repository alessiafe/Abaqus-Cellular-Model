C Engineering constants for tracheids and matrix
      IF (MATNAME == 'MATRIX') THEN
      KW_E1 = (3.407D0)
      KW_E2 = (3.407D0)
      KW_E3 = (3.407D0)
      KW_nu23 = (0.348D0)
      KW_nu13 = (0.348D0)
      KW_nu12 = (0.348D0)
      KW_G23 = (1.266D0)
      KW_G13 = (1.266D0)
      KW_G12 = (1.266D0)
      ELSE IF (MATNAME == 'EW') THEN
      KW_E1 = (49.387D0)
      KW_E2 = (7.887D0)
      KW_E3 = (4.824D0)
      KW_nu23 = (0.353D0)
      KW_nu13 = (0.175D0)
      KW_nu12 = (0.702D0)
      KW_G23 = (1.622D0)
      KW_G13 = (1.649D0)
      KW_G12 = (6.149D0)
      ELSE IF (MATNAME == 'TW') THEN
      KW_E1 = (51.510D0)
      KW_E2 = (6.627D0)
      KW_E3 = (4.804D0)
      KW_nu23 = (0.353D0)
      KW_nu13 = (0.137D0)
      KW_nu12 = (0.807D0)
      KW_G23 = (1.642D0)
      KW_G13 = (1.672D0)
      KW_G12 = (5.979D0)
      ELSE IF (MATNAME == 'LW') THEN
      KW_E1 = (52.764D0)
      KW_E2 = (5.812D0)
      KW_E3 = (4.776D0)
      KW_nu23 = (0.354D0)
      KW_nu13 = (0.104D0)
      KW_nu12 = (0.898D0)
      KW_G23 = (1.655D0)
      KW_G13 = (1.687D0)
      KW_G12 = (5.869D0)
      ELSE
      PRINT *, 'Error: Material ', MATNAME, ' not supported.'
      RETURN
      END IF
