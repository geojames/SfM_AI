# Sample Data

### Point cloud
The required columns/headers are x,y,z
 - Any extra columns like color [r,g,b] will be copied to the output point cloud without any modification.

### Camera positions
The required columns/headers are x,y,z,pitch,roll,yaw
 - Camera/photo names can be included, but are not required
 - This file should be From Agisoft Photoscan/Metashape: Export the Estimated Positions/Orientations from your project
 - Camera Quality can be exported as well

### Sensor properties
The required columns/headers are focal,sensor_x,sensor_y
 - Focal is the focal length of the camera in millimeters
 - sensor_x & sensor_y are the physical sensor dimensions in millimeters
