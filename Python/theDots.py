#!/usr/bin/python2.7	
#
#let the dots (e.g. ants) run around and they'lll do stuff
#
#credit to http://www.petercollingridge.co.uk/pygame-physics-simulation for guidance

import sys
import pygame
import time
import random
import math
from ROOT import TCanvas, TF1, TH1F, TH1I, TFile
import cProfile

#global variables
size_1 = 5 #size of first species
background_colour = (255,255,255) #white
fps = 120



#Basic dot species
class Dot:
  def __init__(self, (x, y),size,init_orientation):
    self.x = x
    self.y = y
    self.size = size
    self.angle = init_orientation
    self.dAngles = [0.,0.,0.,0.,0.] #memory of 5 steps
    self.colour = (255,128,0) #RGB color
    self.thickness = 0 #full dot
  
  def status(self):
    if self.colour == (255,128,0): return 'moving'
    if self.colour == (255,0,0): return 'talking'
    
  def setStatus(self,status):
    if status=='moving': self.colour=(255,128,0)
    if status=='talking': self.colour=(255,0,0)

  def display(self,screen):
    pygame.draw.circle(screen, self.colour, (int(self.x), int(self.y)), self.size, self.thickness)
    
  #movement
  #  change in angle: jitter with "dAngle" around the average delta angle of the last 5 steps
  #  there is a max angle = max turning radius  
  def move(self,screen,step=1,dAngle=0.003,dAnglemax=0.015):
    if self.status()=='talking':
      r = random.random()
      self.x += r*math.sin(self.angle) *step/(0.05*fps) #avoids walking together
      self.y += math.cos(self.angle) *step/(0.05*fps)
      delta_angle = 0
    else:
      average = sum(self.dAngles)/len(self.dAngles)
      delta_angle = random.uniform(average-dAngle,average+dAngle)
      while(abs(delta_angle)>dAnglemax): delta_angle = random.uniform(average-dAngle,average+dAngle)
      self.dAngles.pop(0)
      self.dAngles.append(delta_angle)
      self.angle = self.angle+delta_angle
      self.x += math.sin(self.angle) * step
      self.y += math.cos(self.angle) * step
    self.bounce_wall(screen)
    return (delta_angle)
  
  #deleting individual dots is faster than whiping the whole screen  
  def delete(self,screen):
    pygame.draw.circle(screen, background_colour, (int(self.x), int(self.y)), self.size, self.thickness)
    
  def bounce_wall(self,screen):
    width, height = screen.get_size()
    if self.x > width - self.size:
      self.x = 2*(width - self.size) - self.x
      self.angle = - self.angle

    elif self.x < self.size:
      self.x = 2*self.size - self.x
      self.angle = - self.angle

    if self.y > height - self.size:
      self.y = 2*(height - self.size) - self.y
      self.angle = math.pi - self.angle

    elif self.y < self.size:
      self.y = 2*self.size - self.y
      self.angle = math.pi - self.angle
      
      
# Grid tool
# Only check for possible interactions with dots whcih are in the same sector
# The sector is the key of a dictionary
# Each key holds a list of dots which are present in the class
# this class:
#    * makes the grid
#    * fills the grid
#    * updates the grid after each move

class Grid:
  def __init__(self, screen, (g1,g2)):
    width, height = screen.get_size()
    self.grid={}
    self.size = g1*g2
    self.size_x = g1
    self.size_y = g2
    self.dx=width/(1.*g1)
    self.dy=height/(1.*g2)
    for i in range(self.size_x):
      for j in range(self.size_y):
        self.grid[(i,j)]=[]
  
  #make the grid
  def fill(self,bugs):  #fill the grid
    for bug in bugs:
      key = self.getSector(bug)
      self.grid[key].append(bug)
    print self.grid
      
  def getSector(self,bug):  #get the sector(s) of a bug
    sector_x = int(bug.x/self.dx)
    sector_y = int(bug.y/self.dy)
    return (sector_x, sector_y)
  
#  def getLocalsFromSector(self,bug):  #get all individuals in a certain sector
#    neighbours=[]
#    return neighours
  
#  def getLocals(self,bug): #get all individuals near a specific dot
#    neighbours=[]
#    sectors = getSectors(bug)
#    for sector in sectors
#      neighbours+=getLocalsFromSector[sector]
#    neighbours = list(set(neighbours)) # make unique collection
#    return neighbours
    
#  def updateGrid(self,bug): #update grid, don`t forget to remove a bug if it leaves a sector
#    sectors = getSectors(bug)
#    for sector in sectors
#      self.grid[sector].append(bug)
 #     self.grid[sector]=list(set(self.grid[sector]))
  
# define the Canvas
def MakeWindow(bc,w,h):
  (width, height) = (w, h)
  screen = pygame.display.set_mode((width, height))
  pygame.display.set_caption('the Dots')
  screen.fill(bc)
  return screen
  
#Make a bunch of dots
def MakeDots(number,screen):
  sx, sy = screen.get_size()
  dots=[]
  for i in range(number):
    x = random.randint(size_1+300,sx-size_1-300)
    y = random.randint(size_1+300,sy-size_1-300)
    angle = random.uniform(0,2*math.pi)
    dots.append( Dot((x,y),size_1, angle) )
  return dots
  
  
def check_interact(d1,d2,dic):
  #check if they are close enough to interact
  dx = d1.x - d2.x
  dy = d1.y - d2.y
  distance = math.hypot(dx,dy)
  if distance < (d1.size + d2.size):
    #check in the dictionary if d1 and d2 just met
    if(dic[d1]!=d2 and dic[d2]!=d1 and d1.status()!='talking' and d2.status()!='talking'):
      interact(d1,d2,dic)
      a=2
  else:
    if(dic[d1]==d2 and dic[d2]==d1): #they just interacted, so it's time to get them moving again
      d1.setStatus('moving')
      d2.setStatus('moving')
      
def interact(d1,d2,dic):
  #interaction
  d1.colour=(255,0,0)
  d2.colour=(255,0,0)
  #adjust dictionary
  dic[d1]=d2
  dic[d2]=d1
  #talking status makes the dots slow down
  d1.setStatus('talking')
  d2.setStatus('talking')
  #an interaction makes them change angle
  dAngle=0.3
  d1.angle = random.uniform(d1.angle-dAngle,d1.angle+dAngle)
  d2.angle = random.uniform(d2.angle-dAngle,d2.angle+dAngle)
  
  
    


def main():

  # make my window
  screen_width = 1200
  screen_height = 800
  world = MakeWindow(background_colour,screen_width,screen_height)
  
  #make a bunch of dots
  dots = MakeDots(60,world)
  #dictionary to keep tracks of interactions
  dots_dic = {x:x for x in dots} 
  #make grid which devides the screen and groups the dots
  grid = Grid(world,(10,10))
  grid.fill(dots)
  
 
  #log histograms
  f = TFile("theDotsStats.root","recreate")
  hDAngles = TH1I('hDAngles','hDAngles',200,-0.1,0.1)

  #Run the thing
  clock = pygame.time.Clock()
  running = True
  while running:
    for i, dot in enumerate(dots):
      dot.delete(world)
      value = dot.move(world,0.2) # move a dot
      hDAngles.Fill(value)
      for dot2 in dots[i+1:]: #this loop takes a lot of CPU, should find a way to get interaction candidates more efficient
        check_interact(dot,dot2,dots_dic)
        
      dot.display(world)
      
    pygame.display.update()
    clock.tick(fps) 
       
    for event in pygame.event.get():
      if event.type == pygame.QUIT:
          running = False
          
  
  f.Write()
  f.Close()
  
  
if __name__ == '__main__':
  main()
