from __future__ import annotations
import arcade
try:
    from scipy import constants as scicon
except ModuleNotFoundError:
    G = 6.6743 * 10 ** -11
else:
    G = scicon.G

import math
import random as rand


# Copyright (C) 2023 Bud P. L. Patterson


#This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

#You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 





HEIGHT = 1000
WIDTH = 1000
TITLE = "GRAV"



class Vec2:

    def __init__(self, x: float, y: float):
        """
        Vector type
        direction is in degrees
        :param x:
        :param y:
        """
        self._x = x
        self._y = y
        #print(self)


    @property
    def x(self):
        return self._x

    @property
    def y(self):
        return self._y


    def __getitem__(self, item: int):
        if item == 0:
            return self._x
        elif item == 1:
            return self._y
        else:
            raise IndexError(f"Invalid index '{item}' for Vec2")


    def __add__(self, other: Vec2) -> Vec2:
        if not isinstance(other, Vec2):
            raise TypeError(f"Unsupported operand type(s) for +: 'Vector' and '{type(other)}'")

        return Vec2(self.x + other.x, self.y + other.y)


    def __sub__(self, other: Vec2) -> Vec2:
        if not isinstance(other, Vec2):
            raise TypeError(f"Unsupported operand type(s) for -: 'Vector' and '{type(other)}'")

        return Vec2(self.x - other.x, self.y - other.y)


    def __mul__(self, other: int | float | Vec2) -> Vec2:
        if not isinstance(other, (int, float, Vec2)):
            raise TypeError(f"Unsupported operand type(s) for *: 'Vector' and '{type(other)}'")

        if type(other) != Vec2:
            return Vec2(self.x * other, self.y * other)
        else:
            #We are doing Vec multiplication
            Vec2(self.x * other.x, self.y * other.y)


    def __matmul__(self, other):
        raise NotImplemented


    def dot(self, other: Vec2) -> float:
        if not isinstance(other, (Vec2)):
            raise TypeError(f"Unsupported operand type(s) for 'dot': 'Vector' and '{type(other)}'")

        out = self * other
        out = out.x + out.y
        return out


    def __neg__(self) -> Vec2:
        return Vec2(-self.x, -self.y)


    def __truediv__(self, other: int | float) -> Vec2:
        if not isinstance(other, (int, float)):
            raise TypeError(f"Unsupported operand type(s) for /: 'Vector' and '{type(other)}'")

        return Vec2(self.x / other, self.y / other)


    def __radd__(self, other: Vec2) -> Vec2:
        return self.__add__(other)

    def __rsub__(self, other: Vec2) -> Vec2:
        return self.__sub__(other)

    def __rmul__(self, other: int | float) -> Vec2:
        return self.__mul__(other)

    def __rdiv__(self, other: int | float) -> Vec2:
        return self.__truediv__(other)

    def __abs__(self):
        return math.sqrt(self.x ** 2 + self.y ** 2)

    def __len__(self):
        return 2


    def __str__(self):
        return f"""
Vector M[{self.x}]
       D[{self.y}]
"""





class Phys:

    @staticmethod
    def law_grav(m1: int, m2: int, dist: float, soft:float=0) -> float:
        if dist == 0:
            return 0
        return (G * m1 * m2) / (dist**2 + soft**2)


    @staticmethod
    def law_grav_vec(m1: int, m2: int, pos1: Vec2, pos2: Vec2, soft:float=0) -> Vec2:
        # http://www.scholarpedia.org/article/N-body_simulations_(gravitational)

        out = G * (
            m1 * m2 * (pos1 - pos2)
            /
            ( abs(pos1 - pos2) ** 2 + soft ** 2) ** (3/2)
        )
        return -out


    @staticmethod
    def law_grav_vec_ent(ent1: Entity | arcade.Sprite, ent2: Entity | arcade.Sprite, soft:float=0):
        return Phys.law_grav_vec(ent1.mass, ent2.mass, ent1.pos, ent2.pos, soft)


    @staticmethod
    def law_second(mass: int, accel: Vec2) -> Vec2:
        return mass * accel


    @staticmethod
    def law_accl_from_force(mass: int, force: Vec2) -> Vec2:
        return force / mass




class Entity(arcade.SpriteCircle):

    def __init__(self, mass, color, x: int, y: int, radius=5):
        super().__init__(radius, color, False,)

        self.mass = mass

        #self.change_x
        #self.change_y
        #self.change_angle

        self.center_x = x
        self.center_y = y

        self.vel = Vec2(0, 0)

        self.age = 0

        self.draw_hit_box((0, 255, 0), 1)

        #Have extra properties
        #self.closest =
        self.closest = [None, None]



    @property
    def pos(self) -> Vec2:
        return Vec2(self.center_x, self.center_y)


    @pos.setter
    def pos(self, val: Vec2) -> None:
        self.center_x = val[0]
        self.center_y = val[1]


    def moveByVec2(self, vec: Vec2) -> None:
        self.center_x = self.center_x + math.cos(vec.y) * vec.x
        self.center_y = self.center_y + math.sin(vec.y) * vec.x


    def distance(self, ent: Entity) -> float:
        return math.sqrt(
            (ent.center_x - self.center_x) ** 2 + (ent.center_y - self.center_y) ** 2
        )


    def getForce(self, soft:float=0) -> Vec2:
        global sim

        mainforce: Vec2 = Vec2(0, 0)
        for ent in sim.objects:
            if ent == self: continue
            mainforce += Phys.law_grav_vec_ent(self, ent, soft)

        return mainforce





    def on_update(self, delta_time: float = 1 / 60):
        global sim

        # https://www.cs.princeton.edu/courses/archive/spr01/cs126/assignments/nbody.html
        # http://www.scholarpedia.org/article/N-body_simulations_(gravitational)


        collide = self.collides_with_list(sim.objects)
        if len(collide) > 0:
            for _ in collide:
                self.vel = -self.vel
                _.vel = -_.vel
                # self.mass += _.mass
                # self.width += _.width
                # self.height += _.height
                # _.kill()


        #find closest ent
        #Ent: Entity, Dist: float

        for ent in sim.objects:
            if ent == self:
                continue

            dist = math.sqrt( (ent.center_x - self.center_x) ** 2 + (ent.center_y - self.center_y) ** 2 )
            if self.closest[1] is None:
                self.closest[1] = dist
            else:
                if dist < self.closest[1]:
                    self.closest[0] = ent
                    self.closest[1] = dist
                else:
                    continue


        #Bounce off edge
        if (self.center_x > WIDTH) or (self.center_x < 0) or (self.center_y > HEIGHT) or (self.center_y < 0):
            self.vel = -self.vel


        #Delete ent if too far
        deldist = 10
        if (self.center_x > WIDTH + deldist) or (self.center_x < 0 - deldist) or (self.center_y > HEIGHT + deldist) or (self.center_y < 0 - deldist):
            #self.kill()
            self.vel = Vec2(0, 0)
            self.pos = Vec2(WIDTH / 2, HEIGHT / 2)





        # Soft=10 is ok
        force = self.getForce(20)

        self.vel += Phys.law_accl_from_force(self.mass, force) * delta_time
        self.pos += self.vel * delta_time


        self.age += delta_time




class Sim(arcade.Window):

    def __init__(self, width, height, title):
        super().__init__(width, height, title)
        arcade.set_background_color(arcade.color.BLACK)

        self.objects = arcade.SpriteList()



    def setup(self):
        pass


    def on_draw(self):
        self.clear()

        self.objects.draw()


    def on_update(self, delta_time: float):
        self.objects.on_update(delta_time)


    def on_mouse_press(self, x: int, y: int, button: int, modifiers: int):
        #TODO make entity on click
        #50000
        #themass = rand.randint(25000, 100_000)
        #themass = 25_000_000
        themass = 25_000_000_000_000_000
        #print(themass)
        self.objects.append(Entity(themass, (0,0,255), x, y))
        pass


    def on_key_press(self, symbol: int, modifiers: int):

        if symbol == arcade.key.R:
            self.objects.clear()

        if symbol == arcade.key.D:
            for _ in range(5):
                self.objects.append(Entity(25_000_000_000_000_000, (0,0,255), rand.randint(0, WIDTH), rand.randint(0, HEIGHT)))

        if symbol == arcade.key.L:
            temp = arcade.SpriteList()
            offset_x = 50
            offset_y = 50
            for x in range(WIDTH):
                if x % 100 != 0:
                    continue
                for y in range(HEIGHT):
                    if y % 100 != 0:
                        continue
                    temp.append(Entity(25_000_000_000, (0,255,0), x + offset_x, y + offset_y, ))
            self.objects.extend(temp)



sim = Sim(WIDTH, HEIGHT, TITLE)

def main():

    sim.setup()
    arcade.run()


if __name__ == "__main__":
    main()












