--
-- Author : Jianshu Wu
-- Date: 2024-02-22
--
-- In this shock-fitted solution, we'll take the edge
-- of the grid as the shock location.

print("Begin stagnation line properties extraction.")

jobName = "nonaka"

Db = 14.0e-3
R = Db/2.0

-- Pick up flow solution at final time
fsol = FlowSolution:new{jobName=jobName, dir=".", tindx="last", nBlocks=4} 

f = io.open("stagnation-line-data.dat", 'w')
f:write("#  x                   y                    transrotational temperature                   vibrational temperature   \n") -- Adjust headers

for 

-- Assuming there's a single stagnation point, modify if not
stagPoint = fsol:get_stagnation_point() 
x = stagPoint.x
y = stagPoint.y
pressure = fsol:get_stagnation_pressure() 
temperature = fsol:get_stagnation_temperature()
-- ... Obtain other properties as needed ...

f:write(string.format("%20.12e %20.12e %20.12e %20.12e ...\n", x, y, pressure, temperature, ...)) 

f:close()
print("Done.")