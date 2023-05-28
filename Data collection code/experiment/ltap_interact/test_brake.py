import numpy as np
import matplotlib.pyplot as plt
import math

pedal_position = np.linspace(-1, 0, 51)
# while True:
#     pygame.event.pump()
#     print(self.joystick.get_button(0))
#     print(self.joystick.get_axis(0))
#     time.sleep(1)

brakeCmd = 1.6 + (2.05 * math.log10(
            -1 * self.joystick.get_axis(2) + 1.4) - 1.2) / 0.92
brakeCmd = 1.6 + (2.05 * np.log10(-1 * pedal_position + 1.4) - 1.2) / 0.92

plt.plot(pedal_position, brakeCmd)