# 3D face reconstruction
The 3D geometry of an object can be defined in a computer model by describing the 3D surface of that object. A possible representation of that surface is by means of a 3D point cloud. That is a collection of 3D points that are represented by their 3D coordinates in a given coordinate system.

However, such a 3D point cloud is limited since it doesn’t provide the connectedness between the points. A richer representation is a so-called 3D surface mesh. This is a collection of 3D points that form the corners of planar polygons. The 3D points are called vertices. The polygons are the so-called faces. Often, triangles are used for the polygon since a triangle is the only planar polygon that has no constraints on the 3D positions of the corners. This is in contrast with, for instance, a quadrangle.

The usual way to represent a 3D surface mesh is by enumerating the list of vertices, and then by defining each face by tabulating the three numbers of its vertices.

Stereo vision uses two images of the same object to find the 3D positions of points on the surface of that object. This can provide a partial surface mesh of the object. It is partial since only surface patches that are visible in both images can be reconstructed. To have a more complete surface mesh, the results of other pairs of images, taken from other points of view, must be merged.

Three subjects are imaged from three different points of view. These three images are indicated by ‘left’, ‘middle’, and ‘right’. Each subject is imaged with five different facial expressions. The adjoining figure provides an example. Also given are two sets of calibration images. A checker board has been imaged with different poses. The squares on the checker board has a size of 10x10 mm^2.

Assignment: Design a method for the creation of the 3D surface meshes of the faces of subjects. This is useful, for instance, for face recognition, but also for medical applications. An example is the grading of facial paralysis due to some neural disorder.

Subtasks:
1. The intrinsic and the external parameters of the camera should be obtained from camera calibration for which the checker board images are available.
2. Non-linear lens deformation should be compensated for, and stereo rectification is needed to facilitate the dense stereo matching.
3. Maybe a global colour normalization to make sure that the so-
4. called ‘Constant Brightness Assumption’ holds true. A normalization could be applied with respect to mean and standard deviation of the colour channels.
5. Stereo matching provides so-called disparity maps which is the input for depth map estimation.
6. We don’t want to have the background within the 3D surface mesh. So, you might want to develop and to use a background detection (maybe edge detection in combination with morphological operations).
7. You also may want to create a map which indicates “unreliable disparities”, so that unreliable depth can be removed from the mesh.
8. To merge two 3D surface meshes to one mesh, you may want to use the “ICP algorithm (iterated closest points). Matlab implementations can be found on the internet.
9. A rough assessment of the quality of the mesh is appreciated. How to quantify quality? Discuss your results.

