from __future__ import annotations

import typing

import arcade
import arcade.gui
try:
    from scipy import constants as scicon
except ModuleNotFoundError:
    G = 6.6743 * 10 ** -11
else:
    G = scicon.G

import math
import random as rand
import statistics as stats


# Copyright (C) 2023 Bud P. L. Patterson


#This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

#You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.



HEIGHT = 800
WIDTH = 1000
TITLE = "GRAV"



#TODO
#   Collision !!
#   CLEAN UP CODE, COMMENT, DE-SPEGETIFY
#   ~~Camera, Zoom, and Scrolling
#   ~~Have Camera follow selected object
#   Controllable Entity
#   Mark Units
#   Impliment GUI
#   Look into using GPU for Physical calculations
#   Inelastic Collisions and Heat lost / Energy; Want to conserve energy in system, Heat or something
#   Show total energy in system






Scalar = int | float


class Vec2:

    def __init__(self, x: Scalar = 0, y: Scalar = 0, polar: bool = False):
        """
        2d Vector Type
        :param x:
        :param y:
        """
        assert x is not None, ValueError(f"Vector 'x' cannot be None")
        assert y is not None, ValueError(f"Vector 'y' cannot be None")
        self._x = x  # R when set to polar
        self._y = y  # Theta when set to polar
        self._polar = polar  # Are we using this vector to represent polar co-ordents?
        #print(self)


    @property
    def x(self):
        if not self._polar:
            return self._x
        else:
            return self._x * math.cos(math.radians(self._y))


    @property
    def y(self):
        if not self._polar:
            return self._y
        else:
            return self._x * math.sin(math.radians(self._y))


    def __getitem__(self, item: int):
        if item == 0:
            return self._x
        elif item == 1:
            return self._y
        else:
            raise IndexError(f"Invalid index '{item}' for Vec2")


    def __add__(self, other: Vec2) -> Vec2:
        if not isinstance(other, Vec2):
            raise NotImplementedError(f"Unsupported operand type(s) for +: 'Vector' and '{type(other)}'")

        return Vec2(self.x + other.x, self.y + other.y)


    def __sub__(self, other: Vec2) -> Vec2:
        if not isinstance(other, Vec2):
            raise NotImplementedError(f"Unsupported operand type(s) for -: 'Vector' and '{type(other)}'")

        return Vec2(self.x - other.x, self.y - other.y)


    def __mul__(self, other: int | float | Vec2) -> Vec2:
        #TODO:
        #   Complex Multiplication
        #   https://www.physicsforums.com/threads/multiplying-a-vector-by-a-complex-number.897606/
        #   https://math.stackexchange.com/questions/4060569/intuition-for-multiplying-vectors-by-complex-scalars
        #
        if not isinstance(other, (int, float, Vec2)):
            raise NotImplementedError(f"Unsupported operand type(s) for *: 'Vector' and '{type(other)}'")

        if type(other) != Vec2:
            return Vec2(self.x * other, self.y * other)
        else:
            #We are doing Vec multiplication
            return Vec2(self.x * other.x, self.y * other.y)


    def __matmul__(self, other):
        raise Exception("TODO")


    def __pow__(self, power, modulo=None) -> Vec2 | float | int:
        # FIXME; INNCORRECT
        # https://www.euclideanspace.com/maths/algebra/vectors/vecAlgebra/powers/index.htm
        raise NotImplementedError
        # init = Vec2(self.x, self.y)
        # for times in range(power):
        #     init *= Vec2(self.x, self.y)
        # return init


    def dot(self, other: Vec2) -> Scalar:
        if not isinstance(other, (Vec2)):
            raise NotImplementedError(f"Unsupported operand type(s) for 'dot': 'Vector' and '{type(other)}'")

        out = self * other
        out = out.x + out.y
        return out


    def cross_faux3d(self, other: Vec2) -> Scalar:
        # https://en.wikipedia.org/wiki/Cross_product#Alternative_ways_to_compute
        # https://stackoverflow.com/questions/243945/calculating-a-2d-vectors-cross-product
        # https://www.nagwa.com/en/explainers/175169159270/
        if not isinstance(other, (Vec2)):
            raise NotImplementedError(f"Unsupported operand type(s) for 'cross': 'Vector' and '{type(other)}'")

        return (self.x * other.y) - (self.y * other.x)



    def cross_faux2d(self, other: Vec2) -> Vec2:
        # My own algorithm, taken from its use in `on_draw` in Sim,
        # Have no idea if it is correct or not, and dont have enough of a math background to know, but seems to do the correctish thing
        if not isinstance(other, (Vec2)):
            raise NotImplementedError(f"Unsupported operand type(s) for 'cross': 'Vector' and '{type(other)}'")

        try:
            return ((self / self.mag) + (other / other.mag)) * ((self.mag + other.mag) / 2)
        except ZeroDivisionError:
            #Don't know if I should return (Empty Vec2) when either mag is 0, but rather not have the program crash
            return Vec2()



    def angle_between(self, other: Vec2) -> Scalar:
        """Return the angle between two vectors in radians"""
        # TODO, DO THIS IN POLAR SPACE?
        #FIXME
        if not isinstance(other, (Vec2)):
            raise NotImplementedError(f"Unsupported operand type(s) for 'angle_between': 'Vector' and '{type(other)}'")

        dot = self.dot(other)
        print("A_B_V", dot / self.mag * other.mag, math.radians(dot / self.mag * other.mag))
        return math.acos(math.radians(dot / self.mag * other.mag))


    @property
    def mag(self) -> Scalar:
        if not self._polar:
            return math.sqrt(self.x**2 + self.y**2)
        else:
            return self._x


    @property
    def heading(self) -> Scalar:
        if not self._polar:
            return math.atan2(self.y, self.x)
        else:
            return self._y


    def __neg__(self) -> Vec2:
        return Vec2(-self.x, -self.y)


    def __truediv__(self, other: Scalar) -> Vec2:
        if not isinstance(other, (Scalar)):
            raise NotImplementedError(f"Unsupported operand type(s) for /: 'Vector' and '{type(other)}'")

        return Vec2(self.x / other, self.y / other)


    def __radd__(self, other: Vec2) -> Vec2:
        return self.__add__(other)

    def __rsub__(self, other: Vec2) -> Vec2:
        return self.__sub__(other)

    def __rmul__(self, other: int | float | Vec2) -> Vec2:
        return self.__mul__(other)

    def __rdiv__(self, other: int | float) -> Vec2:
        return self.__truediv__(other)

    def __rdivmod__(self, other: Scalar) -> Vec2:
        return self.__truediv__(other)

    def __abs__(self):
        return math.sqrt(self.x ** 2 + self.y ** 2)

    def __len__(self):
        return 2

    def __str__(self):
        return f"""
Vector <{self._x}, {self._y}> {"P" if self._polar else ""}
"""

    def __repr__(self):
        return str(self)


class MathUtil:

    # Vector Cross product of 2d?
    # Skoggy says (x1 * x2) + (y1 * y2)

    @staticmethod
    def cartesian_to_polar(x: Scalar, y: Scalar) -> tuple[Scalar, Scalar]:
        # https://www.varsitytutors.com/hotmath/hotmath_help/topics/polar-coordinates
        """
        Converts X, Y in cartesian plane to polar co-ords, returns theta in degrees, returns (r, theta)
        :param x: X
        :param y: Y
        :return: (R, θ)
        """
        r = math.sqrt(x ** 2 + y ** 2)
        theta = math.degrees(math.atan(y / x))
        return r, theta


    @staticmethod
    def polar_to_cartesian(r: Scalar, theta: Scalar) -> tuple[Scalar, Scalar]:
        """
        Converts polar R, θ to X, Y cartesian;
        :param r: R
        :param theta: θ, In Degrees
        :return: (X, Y)
        """
        x = r * math.cos(math.radians(theta))
        y = r * math.sin(math.radians(theta))
        return x, y




class Phys:

    @staticmethod
    def law_grav(m1: int, m2: int, dist: float, soft:float=0) -> float:
        if dist == 0:
            return 0
        return (G * m1 * m2) / (dist**2 + soft**2)


    @staticmethod
    def law_grav_vec(m1: int, m2: int, pos1: Vec2, pos2: Vec2, soft:int=0) -> Vec2:
        # http://www.scholarpedia.org/article/N-body_simulations_(gravitational)

        out = G * (
            m1 * m2 * (pos1 - pos2)
            /
            (abs(pos1 - pos2) ** 2 + soft ** 2) ** (3/2)
        )
        return -out


    @staticmethod
    def standard_grav_parameter(m1: int | float, m2: int | float) -> int | float:
        #FIXME, Negative?
        return G * (m1 + m2)


    @staticmethod
    def law_grav_vec_ent(ent1: Entity | arcade.Sprite, ent2: Entity | arcade.Sprite, soft:int=0):
        return Phys.law_grav_vec(ent1.mass, ent2.mass, ent1.pos, ent2.pos, soft)


    @staticmethod
    def law_second(mass: Scalar, accel: Vec2) -> Vec2:
        return mass * accel


    @staticmethod
    def law_accl_from_force(mass: Scalar, force: Vec2) -> Vec2:
        return force / mass


    @staticmethod
    def calc_momentum(mass: Scalar, vel: Vec2) -> Vec2:
        return mass * vel

    @staticmethod
    def elastic_collision(m1: Scalar, m2: Scalar, vel1: Vec2, vel2: Vec2) -> tuple[Vec2, Vec2]:
        """
        Calculates the elastic collision of two entitys and returns the result in a tuple;
        First value is self's new velocity, second is the affected ent
        """
        # https://en.wikipedia.org/wiki/Elastic_collision#One-dimensional_Newtonian

        #TODO: Probably should be using this:
        # https://en.wikipedia.org/wiki/Elastic_collision#Two-dimensional



        #combined mass, because it shows up alot

        smass = abs(m1)
        emass = abs(m2)


        cm = smass + emass

        v1 = \
        (((smass - emass) / cm) * vel1) + (((2 * emass) / cm) * vel2)

        v2 = \
        ((2 * smass) / cm * vel1) + (((emass - smass) / cm) * vel2)

        return (v1, v2)


    @staticmethod
    def inelastic_collision(C: Scalar, m1: Scalar, m2: Scalar, vel1: Vec2, vel2: Vec2) -> tuple[Vec2, Vec2]:
        """
        Calculates the inelastic collision of two entitys and returns the result in a tuple;
        First value is self's new velocity, second is the affected ent

        :param C: This is the Coefficient of Restitution: float 0-1; 1 being completely elastic collision, 0 being completly inelastic collision
        """
        # https://en.wikipedia.org/wiki/Inelastic_collision#Formula


        m1 = abs(m1)
        m2 = abs(m2)

        v1 = \
        (C*m2 * (vel2 - vel1) + (m1*vel1) + (m2*vel2)) / (m1 + m2)

        v2 = \
        (C*m1 * (vel1 - vel2) + (m1*vel1) + (m2*vel2)) / (m1 + m2)

        return (v1, v2)


    @staticmethod
    def elastic_collision_2d():
        # TODO
        """
        Calculates the elastic collision of two entitys and returns the result in a tuple;
        First value is self's new velocity, second is the affected ent
        """
        # https://en.wikipedia.org/wiki/Elastic_collision#Two-dimensional


    @staticmethod
    def calc_orbit_speed_circular(m1: Scalar, m2: Scalar, r: Scalar) -> Scalar:
        """Circular Orbit Speed Calculation"""
        # FIXME POTENTIALLY WRONG
        return math.sqrt(
            (G * (m1 + m2)) / r
        )
        # https://physics.stackexchange.com/questions/546229/what-determines-if-an-object-will-stay-in-a-planets-orbit


    @staticmethod
    def grav_potential_energy(m1: Scalar, m2: Scalar, dist: Scalar) -> Scalar:
        # https://en.wikipedia.org/wiki/Potential_energy#Potential_energy_for_gravitational_forces_between_two_bodies
        # Probably Wrong Tbh
        return - (G * m1 * m2) / dist


    @staticmethod
    def reduced_mass(m1: Scalar, m2: Scalar) -> Scalar:
        # https://en.wikipedia.org/wiki/Reduced_mass
        return m1 * m2 / m1 + m2


    @staticmethod
    def effective_potential():
        # https://en.wikipedia.org/wiki/Effective_potential
        # https://physics.stackexchange.com/questions/83602/determine-if-an-object-is-in-orbit
        pass


    @staticmethod
    def calc_angular_momentum(m: Scalar, r: Vec2, v: Vec2) -> Vec2:
        """

        :param m: Mass of orbiting Body
        :param r: Position Vector
        :param v: Velocity Vector
        :return:
        """
        # https://en.wikipedia.org/wiki/Angular_momentum#Examples
        # http://www.scholarpedia.org/article/Celestial_mechanics#Newton.E2.80.99s_Celestial_Mechanics
        # https://en.wikipedia.org/wiki/Specific_angular_momentum

        #return m * (r.cross_faux2d(v))
        return r.cross_faux2d(m * v)


    @staticmethod
    def calc_specific_angular_momentum(r: Vec2, v: Vec2) -> Vec2:
        """

        :param r: Position Vector
        :param v: Velocity Vector
        :return:
        """
        return r.cross_faux2d(v)


    @staticmethod
    def calc_specific_orbital_energy(m1: Scalar, m2: Scalar, r: Vec2, v: Vec2) -> Scalar:
        # https://en.wikipedia.org/wiki/Specific_orbital_energy
        # https://space.stackexchange.com/questions/1904/how-to-programmatically-calculate-orbital-elements-using-position-velocity-vecto

        # https://space.stackexchange.com/questions/57909/how-can-i-find-the-x-y-position-of-apogee-and-perigee-if-i-only-have-the-value
        """

        :param m1: Mass of object 1
        :param m2: Mass of object 2
        :param r:  Radial Distance Vector, using a relative polar vector
        :param v:  Velocity Vector
        :return:
        """

        return (v.mag ** 2 / 2) - (Phys.standard_grav_parameter(m1, m2) / r.mag)




class Astronomy:
    """
    Class to hold misc stuff, e.g. (Planet names and other stuff)
    """

    # TODO Break down by type?
    planet_names = [
        "Mercury",
        "Venus",
        "Earth",
        "Mars",
        "Jupiter"
    ]


    # https://space.stackexchange.com/questions/1904/how-to-programmatically-calculate-orbital-elements-using-position-velocity-vecto

    # https://space.stackexchange.com/questions/58136/compute-orbital-parameters-from-transverse-and-radial-velocity-at-a-given-distan

    # https://en.wikipedia.org/wiki/Vis-viva_equation

    # https://en.wikipedia.org/wiki/Electrodynamic_tether , not related, just really cool

    @staticmethod
    def calc_orb_period_circle(m1: Scalar, m2: Scalar, rad: Scalar) -> Scalar:
        # https://en.wikipedia.org/wiki/Circular_orbit#Angular_speed_and_orbital_period
        return 2 * math.pi * math.sqrt(
            rad ** 3 / Phys.standard_grav_parameter(m1, m2)
        )


    @staticmethod
    def calc_semi_major_axis(m1: Scalar, m2: Scalar, r: Vec2, v: Vec2) -> Scalar:
        # https://space.stackexchange.com/questions/57909/how-can-i-find-the-x-y-position-of-apogee-and-perigee-if-i-only-have-the-value
        # http://orbitsimulator.com/formulas/OrbitalElements.html
        """

        :param m1: Mass of object 1
        :param m2: Mass of object 2
        :param r:  Radial Distance Vector, using a relative polar vector
        :param v:  Velocity Vector
        :return:
        """
        #return Phys.standard_grav_parameter(m1, m2) / 2 * Phys.calc_specific_orbital_energy(m1, m2, r, v)

        # Formula taken from line 579, of OrbitalElements.html,
        # http://orbitsimulator.com/formulas/OrbitalElements.html ^
        return 1 / (2 / r.mag - v.mag **2 / Phys.standard_grav_parameter(m1, m2))


    @staticmethod
    def calc_eccentricity_vector(m1: Scalar, m2: Scalar, r: Vec2, v: Vec2, d2: bool = False) -> Scalar:
        # https://space.stackexchange.com/questions/57909/how-can-i-find-the-x-y-position-of-apogee-and-perigee-if-i-only-have-the-value
        # THANK YOU RANDOM INTERNET MATH STRANGERS, YOU ARE HEROS!!!

        # Formula taken from line 589, of OrbitalElements.html,
        # http://orbitsimulator.com/formulas/OrbitalElements.html ^
        # Actual God-send, would not be possibe without
        # Thank you, Tony Dunn
        """"""


        e = math.sqrt(1 - (
            r.cross_faux3d(v) ** 2 / Phys.standard_grav_parameter(m1, m2)
        ) / Astronomy.calc_semi_major_axis(m1, m2, r, v))
        return e




class Orbit:

    def __init__(self):
        # https://en.wikipedia.org/wiki/Orbital_elements
        self._eccentricity = 0
        self._semiMajorAxis = 0
        self._periapsis = 0
        self._apoapsis = 0
        #True Anomaly?

    def anomaly(self):
        pass



class Entity(arcade.Sprite):

    def __init__(self, mass: int | float, x: int, y: int, file: str, scale: float=1, *args, **kwargs):

        super().__init__(filename=file, scale=scale, *args, **kwargs)

        #TODO, Determine Entity's orbital state, e.g. (Orbiting, SubOrbital?, "Static")

        self.name = rand.choice(Astronomy.planet_names)
        self.info_text = arcade.Text(
            f"{self.name}",
            (x + (self.width // 2)) * math.cos(math.radians(45)), (y + (self.width // 2)) * math.sin(math.radians(45)),
            font_size=10, anchor_x="left")


        self.mass = mass
        self.file = file

        self.center_x = x
        self.center_y = y

        self.vel = Vec2(0, 0)

        self.age = 0

        #self.draw_hit_box((0, 255, 0), 1)

        #Have extra properties

        self.closest: list[Entity | None, int] = [None, 0]

        self.past_positions: list[tuple[float, float]] = []
        self.trail_color = (rand.randint(0, 255), rand.randint(0, 255), rand.randint(0, 255), rand.randint(150, 255),)

        self.static = False




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


    def getForce(self, soft: int = 0) -> Vec2:
        global sim

        mainforce: Vec2 = Vec2(0, 0)
        for ent in sim.objects:
            if ent == self:
                continue
            mainforce += Phys.law_grav_vec_ent(self, ent, soft) * sim.obj_scale

        return mainforce




    def on_update(self, delta_time: float = 1 / 60):
        global sim

        # https://www.cs.princeton.edu/courses/archive/spr01/cs126/assignments/nbody.html

        self.past_positions.append((self.pos.x, self.pos.y))
        # Prevents list from getting too long
        # if len(self.past_positions) > 300:
        #     self.past_positions.pop(0)

        # if self.age % 60 < .5:
        #     #auto clear path list, as to prevent clutter and slowdown
        #     self.past_positions.clear()


        self.info_text.x = (self.center_x + (self.width // 4))# * math.cos(math.radians(45))
        self.info_text.y = (self.center_y + (self.width // 4))# * math.sin(math.radians(45))

        self.age += delta_time


        collide = self.collides_with_list(sim.objects)
        if len(collide) > 0:
            for _ in collide:
                #inelastic collistion, '.8' works well
                newvels = Phys.elastic_collision(self.mass, _.mass, self.vel, _.vel)

                self.vel = newvels[0]
                _.vel = newvels[1]
                #print(self.closest)

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



        if self.static:
            self.info_text.color = arcade.color.DARK_RED
            #self.vel = Vec2()
            return
        else:
            self.info_text.color = arcade.color.WHITE
            # Probably could write this in a far better way, running on > 5 hours sleep... will do later


        # Soft=10 is ok
        # Soft=20 is better
        # Soft=100 is even better?
        #force = self.getForce(50)
        force = self.getForce(5)

        self.vel += Phys.law_accl_from_force(self.mass, force) * delta_time
        self.pos += self.vel * delta_time



    def __str__(self):
        return f"{self.name} {self.mass / 10 ** 24} TKG, ({self.pos.x}, {self.pos.y})"













from pyglet.math import Vec2 as PVec2
from typing import NewType



#TODO READ THIS
#   https://forum.kerbalspaceprogram.com/index.php?/topic/161861-orbital-calculation-headaches/
#   KSP FOURMS COMING IN CLUTCH



class Sim(arcade.Window):

    def __init__(self, width, height, title):
        super().__init__(width, height, title, resizable=True)
        arcade.set_background_color(arcade.color.BLACK)

        self.cam_ui = arcade.Camera(WIDTH, HEIGHT)
        self.cam_sprites = arcade.Camera(WIDTH, HEIGHT)

        self.game_time = 0

        # WIDTH // 2, HEIGHT // 2
        self.view_pos = Vec2(0, 0)
        self.view_speed = .4


        self.objects = arcade.SpriteList()
        self.center_mass = arcade.SpriteCircle(5, (127, 127, 127))

        self.symbolic_middle = arcade.Sprite(None, center_x=0, center_y=0)
        self.symbolic_middle.visible = False


        #self.manager = arcade.gui.UIManager()
        #self.manager.enable()


        #Placement state ,,, state
        self.place_state_text = arcade.Text(
            "Static", self.width - 5, self.height - 5, arcade.color.WHITE, 12, anchor_x="right", anchor_y="top")

        # Stores how objects should be placed:
        #   S: Static, impart no momentum apon placement
        #   O: Orbit, Orbit nearest body

        self.place_states = ["Static", "Orbit"]
        self.place_state = "Static"



        self.obj_scale = 1

       #self.standard_mass = 25_000_000_000_000_000
        self.standard_mass = 1_000_000_000_000_000

        self.bodys = [
            Entity(self.standard_mass // 2, 0, 0, "media/earthlike.png", ),
            Entity(self.standard_mass, 0, 0, "media/earthlike.png", 2),
            Entity(-self.standard_mass, 0, 0, "media/red.png", ),
            #Entity(self.standard_mass * 5, 0, 0, "media/gas_giant_25.png", .4)
            Entity(self.standard_mass * 500, 0, 0, "media/gas_giant_25.png", .4),
        ]

        #Mouse hover shadow
        self.hover_shadow = self.bodys[0]
        self.hover_shadow.color = (127, 127, 127)
        self.place = True

        for _ in self.bodys:
            _.color = (170, 170, 170)

        self.current_body = 0


        self.draw_vecs = {
            "v": False,
            "t": False,
        }


        self.selected_object: list[Entity | arcade.Sprite] = []
        self.tracking_object = True  #TODO

        self.paused = False
        self.paused_text = arcade.Text(
            "PAUSED", self.width // 2, 2, arcade.color.DARK_PASTEL_RED, 14, anchor_x="center"
        )


        self.m_text = arcade.Text("", 0, 0, font_size=10)  #font size 8


        self.selected_object_readout = arcade.Text(
            "", 0, 150, font_size=10, anchor_x="left", anchor_y="top", multiline=True, width=2250  # 255
        )


        self.w_pressed = False
        self.s_pressed = False
        self.a_pressed = False
        self.d_pressed = False
        self.cam_speed = 2




    def setup(self):
        pass



    def on_draw(self):
        # TODO COMMENT THIS MESS
        self.clear()

        self.cam_sprites.use()

        self.objects.draw()
        for _ in self.objects:
            _.info_text.draw()

        self.center_mass.draw()

        if self.draw_vecs['v']:
            for obj in self.objects:
                arcade.draw_line(
                    obj.center_x, obj.center_y,
                    obj.center_x + obj.vel.x, obj.center_y + obj.vel.y,
                    arcade.color.WHITE, 2
                )

            for obx, obj in enumerate(self.selected_object):
                force = obj.getForce(50)
                try:
                    #forward: Vec2 = ((force / force.mag) + (obj.vel / obj.vel.mag)) * ((force.mag + obj.vel.mag) / 2)
                    forward: Vec2 = force.cross_faux2d(obj.vel)
                except ZeroDivisionError:
                    break

                arcade.draw_line(
                    obj.center_x, obj.center_y,
                    obj.center_x + force.x, obj.center_y + force.y,
                    arcade.color.RED
                )

                # arcade.draw_arc_outline(
                #     obj.center_x, obj.center_y,
                #     forward.x / 100, forward.y / 100,
                #     arcade.color.BLUE,
                #     forward.heading,
                #     -forward.heading
                # )

                arcade.draw_line(
                    obj.center_x, obj.center_y,
                    obj.center_x + forward.x, obj.center_y + forward.y,
                    arcade.color.AMARANTH_PURPLE
                )



        if self.draw_vecs['t']:
            for obx, obj in enumerate(self.selected_object):
                arcade.draw_points(
                    obj.past_positions,
                    # arcade.color.AZURE,
                    #arcade.color.AZURE if obx == 0 else (50, 50, (obx / len(self.selected_object)) * 255)
                    arcade.color.AZURE if obx == 0 else obj.trail_color
                    # (0, 127 - (obx % 6), 255 - (obx % 9))
                    # (round(((obx % 5) / 5) * 255) * 10, 255, round(((obx % 3) / 3) * 255))
                )





        if self.selected_object:
            for obx, obj in enumerate(self.selected_object):
                arcade.draw_circle_outline(
                    obj.center_x, obj.center_y,
                    (max(obj.width, obj.height) // 2),
                    #arcade.color.GREEN
                    #(round(((obx % 5) / 5) * 255) * 10, 255, round(((obx % 3) / 3) * 255) )
                    # ((obx // len(self.selected_object)) * 255, 50, 100)
                    arcade.color.GREEN if obx == 0 else obj.trail_color
                )




        # God Help Me
        # if len(self.selected_object) == 2:
        #     arcade.draw_ellipse_outline(
        #         self.selected_object[0].center_x, self.selected_object[0].center_y,
        #
        #     )







        if len(self.selected_object) > 1:
            mass_center_vec = Vec2(0, 0)

            total_mass = sum(map(lambda x: x.mass, self.selected_object))

            mass_center_vec += sum(map(lambda x: x.pos * x.mass, self.selected_object), start=Vec2(0, 0))

            mass_center_vec /= total_mass
            arcade.draw_circle_filled(
                mass_center_vec.x, mass_center_vec.y, 5, arcade.color.GRAY_BLUE
            )



        if (self.place_state == "Orbit") and self.selected_object and self.place:
            curBody = self.hover_shadow

            # math.sqrt(
            #     ((curBody.center_x + self.cam_sprites.position.x) - self.selected_object.center_x) ** 2
            #     + ((curBody.center_y + self.cam_sprites.position.y) - self.selected_object.center_y) ** 2
            # ),

            arcade.draw_circle_outline(
                self.selected_object[0].center_x, self.selected_object[0].center_y,
                math.sqrt(
                    ((curBody.center_x + self.cam_sprites.position.x) - self.selected_object[0].center_x) ** 2
                    + ((curBody.center_y + self.cam_sprites.position.y) - self.selected_object[0].center_y) ** 2
                ),

                # ((self.cam_sprites.scale - 1) * 100)

                arcade.color.GRAY_BLUE
            )



        self.cam_ui.use()

        self.hover_shadow.draw()

        self.m_text.draw()

        self.place_state_text.draw()

        if self.selected_object:
            self.selected_object_readout.draw()

        if self.paused:
            self.paused_text.draw()






    def on_update(self, delta_time: float):
        self.game_time += delta_time

        #every minute, resync the game scale size to a random object to ensure that floating point error
        #   does not drift the real scale and the measured scale too far

        # if (self.game_time % 60) < 2:
        #     if self.objects:
        #         #self.obj_scale = rand.choice(self.objects).scale
        #         self.obj_scale = stats.mode(map(lambda x: x.scale, self.objects))
        #         print(f"Scale Resync!, new scale: {self.obj_scale}")
        #     else:
        #         print("Failed to resync scale!, no objects to sync to!")


        if not self.paused:
            self.objects.on_update(delta_time)


        # Update place state text
        self.place_state_text.text = self.place_state


        # Find center of mass
        # https://www.khanacademy.org/science/physics/linear-momentum/center-of-mass/v/center-of-mass-equation

        if len(self.objects) > 0:
            self.center_mass.visible = True
        else:
            self.center_mass.visible = False

        mass_center_vec = Vec2(0, 0)
        total_mass = 0
        for obj in self.objects:
            mass_center_vec += obj.pos * obj.mass
            total_mass += obj.mass

        if total_mass > 0:
            mass_center_vec /= total_mass
            self.center_mass.center_x = mass_center_vec.x
            self.center_mass.center_y = mass_center_vec.y





        if self.selected_object:
            self.view_pos = self.selected_object[0].pos
        else:
            if self.w_pressed:
                self.view_pos = Vec2(self.view_pos.x, self.view_pos.y + self.cam_speed)
            if self.s_pressed:
                self.view_pos = Vec2(self.view_pos.x, self.view_pos.y - self.cam_speed)
            if self.a_pressed:
                self.view_pos = Vec2(self.view_pos.x - self.cam_speed, self.view_pos.y)
            if self.d_pressed:
                self.view_pos = Vec2(self.view_pos.x + self.cam_speed, self.view_pos.y)


        if self.selected_object:

            #TODO
            # CLEAN UP THIS MESS. HOLY SHIT.
            systemtext = ""
            if len(self.selected_object) > 1:

                needed_transverse_speed = round(Phys.calc_orbit_speed_circular(self.selected_object[0].mass, self.selected_object[1].mass,
                          arcade.get_distance_between_sprites(self.selected_object[0], self.selected_object[1])), 2)

                orbiting = (needed_transverse_speed * .9 <= round(self.selected_object[1].vel.mag) <= needed_transverse_speed * 1.1)
                # https://physics.stackexchange.com/a/546233
                # TODO
                #  CLEAN UP THIS CODE,
                #  MOVE THIS TEXT GENERATION OR ATLEAST THE PARAMETER GENERATION TO A SEPARATE METHOD!!!


                polar_co_ords = MathUtil.cartesian_to_polar(
                    self.selected_object[1].pos.x - self.selected_object[0].pos.x,
                    self.selected_object[1].pos.y - self.selected_object[0].pos.y,
                )
                polar_co_ords_vec2 = Vec2(polar_co_ords[0], polar_co_ords[1], True)
                diff_co_vec2 = Vec2(
                    self.selected_object[1].pos.x - self.selected_object[0].pos.x,
                    self.selected_object[1].pos.y - self.selected_object[0].pos.y
                )

                systemtext = f"""
{' - '.join(map(lambda x: x.name, self.selected_object))} System:
Total Mass: {sum(map(lambda x: x.mass, self.selected_object)) / 10 ** 12:,} TU
Gravitational Potential Energy: {
                round(Phys.grav_potential_energy(self.selected_object[0].mass, self.selected_object[1].mass,
                                                 arcade.get_distance_between_sprites(self.selected_object[0], self.selected_object[1])) / 10 ** 12 * self.obj_scale, 2):,} TJ
Reduced Mass: {Phys.reduced_mass(self.selected_object[0].mass, self.selected_object[1].mass) / 10 ** 12:,} TU
Distance: {round(arcade.get_distance_between_sprites(self.selected_object[0], self.selected_object[1]), 2):,}
SEMI-MAJOR-AXIS: {Astronomy.calc_semi_major_axis(
self.selected_object[1].mass, self.selected_object[0].mass, polar_co_ords_vec2, self.selected_object[1].vel
):,}
SOE: {Phys.calc_specific_orbital_energy(
self.selected_object[1].mass, self.selected_object[0].mass, polar_co_ords_vec2, self.selected_object[1].vel                    
):,}
Angular Momentum: {round(
Phys.calc_angular_momentum(self.selected_object[1].mass, self.selected_object[1].pos, self.selected_object[1].vel).mag / 10 ** 12, 2):,}
Radial Distance Vector: {
polar_co_ords_vec2
}
Ecc: {
Astronomy.calc_eccentricity_vector(self.selected_object[0].mass, self.selected_object[1].mass, polar_co_ords_vec2, self.selected_object[1].vel):,}
Vel (1): {
self.selected_object[1].vel
}
"""


            # Make mass readout in Kilo-Units, Round age?,
            individual_text = f"""
{self.selected_object[0].name}:
Mass: {self.selected_object[0].mass / 10 ** 12 :,} TU
X, Y: {round(self.selected_object[0].center_x, 2)}, {round(self.selected_object[0].center_y, 2)}
Age:  {round(self.selected_object[0].age, 2)}
Velocity: {round(self.selected_object[0].vel.x, 1), round(self.selected_object[0].vel.y, 1)}
Speed: {round(self.selected_object[0].vel.mag)}
"""


            if len(self.selected_object) == 1:
                self.selected_object_readout.text = individual_text
            else:
                self.selected_object_readout.text = systemtext

            self.selected_object_readout.y = self.selected_object_readout.content_height
            #self.selected_object_readout.width = self.selected_object_readout.content_width + 10

        else:
            self.selected_object_readout.text = ""




        #FIXME, Scale the shadow in relation to the sprites cam, might just have to make another cam for it
        #self.hover_shadow.scale = 1 / self.cam_sprites.scale



        self.scroll()


    def scroll(self):
        pos = PVec2(
            (self.view_pos.x * 1 / self.cam_sprites.scale) - self.width / 2,
            (self.view_pos.y * 1 / self.cam_sprites.scale) - self.height / 2
        )
        self.cam_sprites.move_to(pos, self.view_speed)




    def on_resize(self, width: float, height: float):
        self.place_state_text.x, self.place_state_text.y = self.width - 5, self.height - 5
        self.paused_text.x = self.width // 2

        self.cam_sprites.resize(int(width), int(height))
        self.cam_ui.resize(int(width), int(height))

        self.view_pos = Vec2(0, 0)
        super().on_resize(width, height)



    def on_mouse_press(self, x: int, y: int, button: int, modifiers: int):

        # print("MOUSE", x, y, "VIEW_POS", self.view_pos, "DMOUSE", self.view_pos.x - x, self.view_pos.y - y)

        #print(self.cam_sprites.position)

        # Updates the X Y Co-ords to the correct screen position by adjusting based on where the camera is looking
        x, y = \
        (self.cam_sprites.position.x) + x, \
        (self.cam_sprites.position.y) + y



        # print("XY", x, y)


        if (button == arcade.MOUSE_BUTTON_LEFT) and (self.place):
            cbody = self.bodys[self.current_body]
            self.objects.append(Entity(cbody.mass, x, y, cbody.file, cbody.scale * self.obj_scale))

            #If we are set to make this body orbit another, find angle and velocity needed to make that happen
            # TODO, make it orbit nearest body if not selected?
            if (self.place_state == "Orbit") and self.selected_object:

                selObj: Entity | arcade.Sprite = self.selected_object[0]
                curObj = self.objects[-1]
                dist = math.sqrt( ((selObj.center_x - curObj.center_x)**2) + ((selObj.center_y - curObj.center_y)**2))

                speed_needed = Phys.calc_orbit_speed_circular(
                    selObj.mass, curObj.mass,
                    dist
                ) * self.obj_scale

                angle = math.atan2(selObj.pos.y - curObj.pos.y, selObj.pos.x - curObj.pos.x)

                #speed_needed /= 100

                #print("ANGLE", angle, "SPEED", speed_needed)

                # θ = Theta
                # Y = SIN θ
                # X = COS θ

                # Turn angle into actual vector
                y_comp = (dist * math.sin(angle))
                x_comp = (dist * math.cos(angle))

                #print("x", x_comp, "y", y_comp)

                new_vel = Vec2(x_comp, y_comp)
                new_vel = Vec2(new_vel.y * -1, new_vel.x)  # Rotate vector by -90 Degrees
                new_vel = new_vel / new_vel.mag  # Need to make it into a unit vector before applying speed
                new_vel *= speed_needed

                curObj.vel = new_vel




        elif (button == arcade.MOUSE_BUTTON_LEFT) and not self.place:
            #Run selecting object code here
            clicked_object = arcade.get_sprites_at_point((x, y), self.objects)

            # Check to see if we actually clicked on an object: If yes make that the current object, if not, make None
            if clicked_object:
                if not modifiers & arcade.key.MOD_SHIFT:
                    # https://api.arcade.academy/en/latest/keyboard.html#keyboard-modifiers
                    self.selected_object.clear()
                    self.selected_object.insert(0, clicked_object[-1])  #-1 to grab the top-most object
                else:
                    if clicked_object[-1] not in self.selected_object:
                        self.selected_object.append(clicked_object[-1])
            else:
                self.selected_object.clear()
            #print(self.selected_object)



        elif button == arcade.MOUSE_BUTTON_RIGHT:
            self.place = not self.hover_shadow.visible
            self.on_mouse_enter(0, 0)



    def on_mouse_enter(self, x: int, y: int):
        self.hover_shadow.visible = self.place

    def on_mouse_leave(self, x: int, y: int):
        self.hover_shadow.visible = False

    def on_mouse_motion(self, x: int, y: int, dx: int, dy: int):
        self.hover_shadow.center_x = x
        self.hover_shadow.center_y = y

        self.m_text.x = x
        self.m_text.y = y

        #===
        nx, ny = \
        (self.cam_sprites.position.x) + x, \
        (self.cam_sprites.position.y) + y


        self.m_text.text = f"({x}, {y}) | ({round(nx, 2)}, {round(ny, 2)})"
        # ({x * self.cam_sprites.scale}, {y * self.cam_sprites.scale})


    def scale(self, obj: Entity | arcade.Sprite, scaler: int | float):
        # obj.rescale_relative_to_point(
        #     ((self.width // 2) + self.cam_sprites.position.x, (self.height // 2) + self.cam_sprites.position.y),
        #     scaler
        # )
        obj.rescale_relative_to_point(
            (0, 0),
            scaler
        )

    def on_key_press(self, symbol: int, modifiers: int):


        if (symbol == arcade.key.I) and (self.selected_object):
            for _ in self.selected_object:
                print(_.pos)

        incri = .125

        if symbol == arcade.key.UP:
            self.obj_scale *= 1 - incri
            for obj in self.objects:
                self.scale(obj, (1 - incri))

        if symbol == arcade.key.DOWN:
            self.obj_scale *= 1 + incri
            for obj in self.objects:
                self.scale(obj, (1 + incri))


        if (symbol == arcade.key.P) and (modifiers != arcade.key.MOD_SHIFT):
            # Reset the scale
            for obj in self.objects:
                self.scale(obj, 1 / self.obj_scale)
            self.obj_scale = 1
        elif (symbol == arcade.key.P) and (modifiers == arcade.key.MOD_SHIFT):
            for obj in self.objects:
                obj.scale = 1



        if symbol == arcade.key.Z:
            for obj in self.selected_object:
                obj.static = not obj.static

        if symbol == arcade.key.C:
            if not modifiers & arcade.key.MOD_SHIFT:
                for obj in self.selected_object:
                    obj.past_positions.clear()
            else:
                for obj in self.objects:
                    obj.past_positions.clear()



        #Move Cam
        if symbol == arcade.key.W:
            self.w_pressed = True
        if symbol == arcade.key.S:
            self.s_pressed = True
        if symbol == arcade.key.A:
            self.a_pressed = True
        if symbol == arcade.key.D:
            self.d_pressed = True



        #Pause
        if symbol == arcade.key.SPACE:
            self.paused = not self.paused


        #Clear all objects
        if symbol == arcade.key.R:
            self.objects.clear()
            self.obj_scale = 1


        # Deletes selected object(s)
        if (symbol == arcade.key.DELETE) and self.selected_object:
            for obj in self.selected_object:
                obj.kill()
            self.selected_object.clear()


        #Show Vectors
        if symbol == arcade.key.V:
            self.draw_vecs['v'] = not self.draw_vecs['v']

        #Show Trails
        if symbol == arcade.key.T:
            self.draw_vecs['t'] = not self.draw_vecs['t']


        #Spawn random earthlikes
        if symbol == arcade.key.M:
            for _ in range(5):
                self.objects.append(Entity(25_000_000_000_000_000,
                                           rand.randint(0, round(self.width + self.cam_sprites.position.x)),
                                           rand.randint(0, round(self.height + self.cam_sprites.position.y)),
                                           "media/earthlike.png", scale=self.obj_scale))


        #Stop all objects
        if symbol == arcade.key.L:
            for _ in self.objects:
                _.vel = Vec2(0, 0)



        #Change placement state
        if symbol == arcade.key.TAB:
            self.place_state = self.place_states[(self.place_states.index(self.place_state) + 1) % len(self.place_states)]


        # if symbol == arcade.key.UP:
        #
        #     print(self.cam_sprites.scale)
        #
        #     if self.cam_sprites.scale < 2:
        #         self.cam_sprites.scale = round(self.cam_sprites.scale + .1, 1)
        #         self.cam_sprites.set_projection()
        #
        #
        # if symbol == arcade.key.DOWN:
        #
        #     print(self.cam_sprites.scale)
        #
        #     if self.cam_sprites.scale > .1:
        #         self.cam_sprites.scale = round(self.cam_sprites.scale - .1, 1)
        #         self.cam_sprites.set_projection()
        #
        #
        # if symbol == arcade.key.P:
        #     self.cam_sprites.scale = 1.0
        #     self.cam_sprites.set_projection()




        #Change current object
        if symbol == arcade.key.BRACKETRIGHT:
            pos = self.bodys[self.current_body].pos
            vis = self.bodys[self.current_body].visible

            self.current_body = (self.current_body + 1) % len(self.bodys)
            self.hover_shadow = self.bodys[self.current_body]

            self.hover_shadow.visible = vis
            self.hover_shadow.pos = pos

        if symbol == arcade.key.BRACKETLEFT:
            pos = self.bodys[self.current_body].pos
            vis = self.bodys[self.current_body].visible

            self.current_body = (self.current_body - 1) % len(self.bodys)
            self.hover_shadow = self.bodys[self.current_body]

            self.hover_shadow.visible = vis
            self.hover_shadow.pos = pos


    def on_key_release(self, symbol: int, modifiers: int):

        #Move Cam
        if symbol == arcade.key.W:
            self.w_pressed = False
        if symbol == arcade.key.S:
            self.s_pressed = False
        if symbol == arcade.key.A:
            self.a_pressed = False
        if symbol == arcade.key.D:
            self.d_pressed = False



sim = Sim(WIDTH, HEIGHT, TITLE)

def main():

    sim.setup()
    arcade.run()


if __name__ == "__main__":
    main()
