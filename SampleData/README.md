# Sample Data
All file header values should be in *lowercase*


### Point cloud
The required columns/headers are x,y,z
 - Any extra columns like color [r,g,b] will be copied to the output point cloud without any modification.
 - Import LAS to CloudCompare, Export CSV file (headers are fine)

### Camera positions
The required columns/headers are x,y,z,pitch,roll,yaw (editing required)
 - Camera/photo names can be included, but are not required
 - This file should be From Agisoft Photoscan/Metashape: Export the Estimated Positions/Orientations from your project
 - Camera Quality can be exported as well
 1) switch to detail view in the Photo windows in Metashape
 2) select all photos, right-click > "Estimate Quality"
 3) Right-click, Copy
 4) paste to excel, copy quality column to exported cam posiiton file

### Sensor properties
The required columns/headers are focal,sensor_x,sensor_y
 - Focal is the focal length of the camera in millimeters
 - sensor_x & sensor_y are the physical sensor dimensions in millimeters
