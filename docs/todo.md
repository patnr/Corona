Todo
================================================
- Note: In order to sum-to-1: omit Q_hosp
- Note: Bug?

         - I2Q_sevr should use dt_sevr, not dt_hosp 
         - d_hosp should switch dt's

- Note: Questionable definition: 

        def Hospitalized(self):
            return self.Q_hosp + self.Q_fatl


