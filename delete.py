import pygame
pygame.init()

R = 100
w = 10
r = 0
g = 0
b = 255
x = 640
y = 480

screen = pygame.display.set_mode((x, y)) #x and y are height and width

pygame.draw.circle(screen, (r,g,b), (x, y), R, w) #(r, g, b) is color, (x, y) is center, R is radius and w is the thickness of the circle border

for i in range(10):
    i += 1