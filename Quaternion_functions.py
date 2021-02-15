import numpy as np

def euler_to_quaternion(roll, pitch, yaw):
    # This function is used to translate the euler angles (roll, pitch, yaw) to quaternions

    qx = np.sin(roll/2) * np.cos(pitch/2) * np.cos(yaw/2) - np.cos(roll/2) * np.sin(pitch/2) * np.sin(yaw/2)
    qy = np.cos(roll/2) * np.sin(pitch/2) * np.cos(yaw/2) + np.sin(roll/2) * np.cos(pitch/2) * np.sin(yaw/2)
    qz = np.cos(roll/2) * np.cos(pitch/2) * np.sin(yaw/2) - np.sin(roll/2) * np.sin(pitch/2) * np.cos(yaw/2)
    qw = np.cos(roll/2) * np.cos(pitch/2) * np.cos(yaw/2) + np.sin(roll/2) * np.sin(pitch/2) * np.sin(yaw/2)

    return [qx, qy, qz, qw]

def quaternion_error(current_quaternion, command_quaternion):
    # For the control of the ADCS the control system will provide a command quaternion. The difference
    # between the current quaternion and the command quaternion is required to produce a change in the 
    # system.

    qc1, qc2, qc3, qc4 = command_quaternion
    q_c = np.array(([qc4, qc3, -qc2, -qc1],[-qc3, qc4, qc1, -qc2],[qc2, -qc1, qc4, -qc3], [qc1, qc2, qc3, qc4]))
    q_error = np.matmul(q_c, current_quaternion)

    return q_error