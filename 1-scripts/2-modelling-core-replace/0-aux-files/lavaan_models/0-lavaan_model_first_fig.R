


#Non-imaging lavaan model. 

model<- ' # direct effects
             CRP ~ c*CM
             
           # mediator 1
             BMI ~ a*CM
             CRP ~ b*BMI
            
          # mediator 2 
            AT ~ a2*CM
            CRP ~ b2*AT
          
           # indirect effect (a*b)
             indirect1 := a*b
             indirect2 := a2*b2
             
           # total effect
             totalCRP := c + indirect1 + indirect2 

'