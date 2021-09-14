"""
Constants
=========
"""
import  numpy as np

N =  96
A = 8 # depth and height of the bounding box
B = 20  # length of the bounding box


def _make_rotation_matrices():
    """
        prepares rotation matrices

    :return: returns 6 matrices with rotations around X,Y,Z axis on plus/minus 90 degrees.
    :rtype: (3,3,6) numpy array
    """

    rot_z_1 = np.array([0, -1, 0, 1, 0, 0, 0, 0, 1]).reshape(3, 3, 1)
    rot_z_2 = np.array([0, 1, 0, -1, 0, 0, 0, 0, 1]).reshape(3, 3, 1)

    rot_x_1 = np.array([1, 0, 0, 0, 0, -1, 0, 1, 0]).reshape(3, 3, 1)
    rot_x_2 = np.array([1, 0, 0, 0, 0, 1, 0, -1, 0]).reshape(3, 3, 1)

    rot_y_1 = np.array([0, 0, 1, 0, 1, 0, -1, 0, 0]).reshape(3, 3, 1)
    rot_y_2 = np.array([0, 0, -1, 0, 1, 0, 1, 0, 0]).reshape(3, 3, 1)

    rot = np.concatenate([rot_x_1, rot_x_2, rot_y_1, rot_y_2, rot_z_1, rot_z_2], axis=2)

    return rot

rot = _make_rotation_matrices()