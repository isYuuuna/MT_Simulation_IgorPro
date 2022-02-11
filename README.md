# N Parallel MTs Simulation
The repository is for storing one-sided microtubules simulation code written in Igor pro. This computational model will help to understand MT bunldes and their behvaiours. 
Due to Github file-size restriction, only procedure files (.ipf) are up-to-date.
The programming language used in Igor pro is similar to Java. <br/>
<br/>
Feb.10 update: Implement change point detection in Igor Pro based on ED-PLET method. When I was trying to find an algorithm for detecting switching point in kinetochore velocity, I learned about change point detection. As I came across this amazing blog: https://aakinshin.net/posts/edpelt/, I thought I could try the same thing in Igor. This algorithm works really well on kinetochore velocity data. I am currently trying to generalize the algorithm (such as trying different penalty and cost function) in order to make it more flexible under different combination of kappa and external forces. 

