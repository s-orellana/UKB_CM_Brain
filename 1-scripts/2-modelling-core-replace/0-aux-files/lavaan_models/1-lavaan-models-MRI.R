

####################  FULL MODEL 
model.full<-'
              #-----------Left side of the model
              #direct effects
               CRP ~ c*CM
               
              #mediator 1
               BMI ~ a*CM
               CRP ~ b*BMI
              
              #mediator 2 
               AT ~ a2*CM
               CRP ~ b2*AT
               
              
              #--------------Brain side of the model
              #CRP
              grey_matter  ~ f*CRP
              
              #BMI & AT
              grey_matter  ~ e*BMI
              grey_matter  ~ g*AT
              
              #--------------INDIRECT effects 
              
              med1 := a*e
              med2 := a*b*f
              med3 := a2*g
              med4 := a2*b2*f
              med5 := c*f
       
                    '



#---- (brain) side
#a*e     := a*e
#a*b*f   := a*b*f
#a2*g    := a2*g
#a2*b2*f := a2*b2*f
#c*f     := c*f

####################  SPARSE MODEL 
model.sparse<-'
              #-----------Left side of the model
              #direct effects
               CRP ~ c*CM
               
              #mediator 1
               BMI ~ a*CM
               CRP ~ b*BMI
              
              #mediator 2 
               AT ~ a2*CM
               CRP ~ b2*AT
               
              
              #--------------Brain side of the model
              #CRP
              grey_matter ~ f*CRP

              
              #--------------INDIRECT effects 

              
              #---- (brain) side
              med2  := a*b*f
              med4  := a2*b2*f
              med5  := c*f

            '