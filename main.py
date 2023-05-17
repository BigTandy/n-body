from __future__ import annotations
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




class Vec2:

    def __init__(self, x: float, y: float):
        """
        Vector type
        direction is in degrees
        :param x:
        :param y:
        """
        assert x is not None, ValueError(f"Vector 'x' cannot be None")
        assert y is not None, ValueError(f"Vector 'y' cannot be None")
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
            raise NotImplemented(f"Unsupported operand type(s) for +: 'Vector' and '{type(other)}'")

        return Vec2(self.x + other.x, self.y + other.y)


    def __sub__(self, other: Vec2) -> Vec2:
        if not isinstance(other, Vec2):
            raise NotImplemented(f"Unsupported operand type(s) for -: 'Vector' and '{type(other)}'")

        return Vec2(self.x - other.x, self.y - other.y)


    def __mul__(self, other: int | float | Vec2 ) -> Vec2:
        #TODO:
        #   Complex Multiplication
        #   https://www.physicsforums.com/threads/multiplying-a-vector-by-a-complex-number.897606/
        #   https://math.stackexchange.com/questions/4060569/intuition-for-multiplying-vectors-by-complex-scalars
        #   
        if not isinstance(other, (int, float, Vec2)):
            raise NotImplemented(f"Unsupported operand type(s) for *: 'Vector' and '{type(other)}'")

        if type(other) != Vec2:
            return Vec2(self.x * other, self.y * other)
        else:
            #We are doing Vec multiplication
            Vec2(self.x * other.x, self.y * other.y)


    def __matmul__(self, other):
        raise Exception("TODO")


    def dot(self, other: Vec2) -> float:
        if not isinstance(other, (Vec2)):
            raise NotImplemented(f"Unsupported operand type(s) for 'dot': 'Vector' and '{type(other)}'")

        out = self * other
        out = out.x + out.y
        return out


    def mag(self) -> float | int:
        return math.sqrt(self.x**2 + self.y**2)


    def heading(self) -> float | int:
        return math.atan2(self.y, self.x)


    def __neg__(self) -> Vec2:
        return Vec2(-self.x, -self.y)


    def __truediv__(self, other: int | float) -> Vec2:
        if not isinstance(other, (int, float)):
            raise NotImplemented(f"Unsupported operand type(s) for /: 'Vector' and '{type(other)}'")

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
Vector <{self.x}, {self.y}>
"""



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
    def law_grav_vec_ent(ent1: Entity | arcade.Sprite, ent2: Entity | arcade.Sprite, soft:int=0):
        return Phys.law_grav_vec(ent1.mass, ent2.mass, ent1.pos, ent2.pos, soft)


    @staticmethod
    def law_second(mass: int, accel: Vec2) -> Vec2:
        return mass * accel


    @staticmethod
    def law_accl_from_force(mass: int, force: Vec2) -> Vec2:
        return force / mass


    @staticmethod
    def elastic_collision(m1: int, m2: int, vel1: Vec2, vel2: Vec2) -> tuple[Vec2, Vec2]:
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
    def inelastic_collision(C: float, m1: int, m2: int, vel1: Vec2, vel2: Vec2) -> tuple[Vec2, Vec2]:
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
    def calc_needed_orbit_speed(m1: int, m2: int, r: float) -> float:
        return math.sqrt(
            (G * (m1 + m2)) / r
        )


    @staticmethod
    def grav_potential_energy(m1: int | float, m2: int | float, dist: int | float):
        # https://en.wikipedia.org/wiki/Potential_energy#Potential_energy_for_gravitational_forces_between_two_bodies
        return - (G * m1 * m2) / dist


    @staticmethod
    def reduced_mass(m1: int | float, m2: int | float) -> int | float:
        # https://en.wikipedia.org/wiki/Reduced_mass
        return m1 * m2 / m1 + m2


    @staticmethod
    def effective_potential():
        # https://en.wikipedia.org/wiki/Effective_potential
        # https://physics.stackexchange.com/questions/83602/determine-if-an-object-is-in-orbit
        pass




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





class Entity(arcade.Sprite):

    def __init__(self, mass: int, x: int, y: int, file: str, scale: float=1, *args, **kwargs):

        super().__init__(filename=file, scale=scale, *args, **kwargs)

        #TODO, Determine Entity's orbital state, e.g. (Orbiting, SubOrbital?, "Static")

        self.name = rand.choice(Astronomy.planet_names)
        self.info_text = arcade.Text(
            f"{self.name}",
            (x + (self.width // 2)) * math.cos(math.radians(45)), (y + (self.width // 2)) * math.sin(math.radians(45)),
            font_size=8, anchor_x="left")


        self.mass = mass
        self.file = file

        self.center_x = x
        self.center_y = y

        self.vel = Vec2(0, 0)

        self.age = 0

        #self.draw_hit_box((0, 255, 0), 1)

        #Have extra properties

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


    def getForce(self, soft: int = 0) -> Vec2:
        global sim

        mainforce: Vec2 = Vec2(0, 0)
        for ent in sim.objects:
            if ent == self:
                continue
            mainforce += Phys.law_grav_vec_ent(self, ent, soft)

        return mainforce




    def on_update(self, delta_time: float = 1 / 60):
        global sim

        # https://www.cs.princeton.edu/courses/archive/spr01/cs126/assignments/nbody.html


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




        # Soft=10 is ok
        # Soft=20 is better
        # Soft=100 is even better?
        force = self.getForce(50)

        self.vel += Phys.law_accl_from_force(self.mass, force) * delta_time
        self.pos += self.vel * delta_time

        self.age += delta_time

        #Update Flavor? Text POS
        self.info_text.x = (self.center_x + (self.width // 4)) * math.cos(math.radians(45))
        self.info_text.y = (self.center_y + (self.width // 4)) * math.sin(math.radians(45))






from pyglet.math import Vec2 as PVec2


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

        #Mouse hover shadow
        self.hover_shadow = Entity(0, 0, 0, "media/earthlike.png", )
        self.hover_shadow.color = (127, 127, 127)
        self.place = True

        self.obj_scale = 1

       #self.standard_mass = 25_000_000_000_000_000
        self.standard_mass = 1_000_000_000_000_000

        self.bodys = [
            Entity(self.standard_mass // 2, 0, 0, "media/earthlike.png", ),
            Entity(-self.standard_mass, 0, 0, "media/red.png", ),
            #Entity(self.standard_mass * 5, 0, 0, "media/gas_giant_25.png", .4)
            Entity(self.standard_mass * 500, 0, 0, "media/gas_giant_25.png", .4)
        ]

        for _ in self.bodys:
            _.color = (170, 170, 170)

        self.current_body = 0


        self.draw_vecs = False


        self.selected_object: Entity | None = None
        self.tracking_object = True  #TODO

        self.paused = False
        self.paused_text = arcade.Text(
            "PAUSED", self.width // 2, 2, arcade.color.DARK_PASTEL_RED, 14, anchor_x="center"
        )


        self.m_text = arcade.Text("", 0, 0, font_size=10)  #font size 8


        self.selected_object_readout = arcade.Text(
            "", 0, 150, font_size=10, anchor_x="left", anchor_y="top", multiline=True, width=225
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
        # for _ in self.objects:
        #     _.info_text.draw()

        self.center_mass.draw()

        if self.draw_vecs:
            for _ in self.objects:
                arcade.draw_line(
                    _.center_x, _.center_y,
                    _.center_x + _.vel.x, _.center_y + _.vel.y,
                    arcade.color.WHITE, 2
                )


        if self.selected_object:
            arcade.draw_circle_outline(
                self.selected_object.center_x, self.selected_object.center_y,
                (max(self.selected_object.width, self.selected_object.height) // 2),
                arcade.color.GREEN
            )



        if (self.place_state == "Orbit") and (self.selected_object is not None) and self.place:
            curBody = self.hover_shadow

            # math.sqrt(
            #     ((curBody.center_x + self.cam_sprites.position.x) - self.selected_object.center_x) ** 2
            #     + ((curBody.center_y + self.cam_sprites.position.y) - self.selected_object.center_y) ** 2
            # ),

            arcade.draw_circle_outline(
                self.selected_object.center_x, self.selected_object.center_y,
                math.sqrt(
                    ((curBody.center_x + self.cam_sprites.position.x) - self.selected_object.center_x) ** 2
                    + ((curBody.center_y + self.cam_sprites.position.y) - self.selected_object.center_y) ** 2
                ),

                # ((self.cam_sprites.scale - 1) * 100)

                arcade.color.GRAY_BLUE
            )



        self.cam_ui.use()

        self.hover_shadow.draw()

        self.m_text.draw()

        self.place_state_text.draw()

        if self.selected_object is not None:
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





        if self.selected_object is not None:
            self.view_pos = self.selected_object.pos
        else:
            if self.w_pressed:
                self.view_pos = Vec2(self.view_pos.x, self.view_pos.y + self.cam_speed)
            if self.s_pressed:
                self.view_pos = Vec2(self.view_pos.x, self.view_pos.y - self.cam_speed)
            if self.a_pressed:
                self.view_pos = Vec2(self.view_pos.x - self.cam_speed, self.view_pos.y)
            if self.d_pressed:
                self.view_pos = Vec2(self.view_pos.x + self.cam_speed, self.view_pos.y)


        if self.selected_object is not None:
            # Make mass readout in Kilo-Units, Round age?,
            self.selected_object_readout.text = f"""
Name: {self.selected_object.name}
Mass: {self.selected_object.mass / 10 ** 12 :,} TU
X, Y: {round(self.selected_object.center_x, 2)}, {round(self.selected_object.center_y, 2)}
Age:  {round(self.selected_object.age, 2)}
Velocity: {round(self.selected_object.vel.x, 1), round(self.selected_object.vel.y, 1)}"""
            self.selected_object_readout.y = self.selected_object_readout.content_height
            #self.selected_object_readout.width = self.selected_object_readout.content_width + 10


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
            if (self.place_state == "Orbit") and (self.selected_object is not None):

                selObj: Entity | arcade.Sprite = self.selected_object
                curObj = self.objects[-1]
                dist = math.sqrt( ((selObj.center_x - curObj.center_x)**2) + ((selObj.center_y - curObj.center_y)**2))

                speed_needed = Phys.calc_needed_orbit_speed(
                    selObj.mass, curObj.mass,
                    dist
                )

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
                new_vel = new_vel / new_vel.mag()  # Need to make it into a unit vector before applying speed
                new_vel *= speed_needed

                curObj.vel = new_vel




        elif (button == arcade.MOUSE_BUTTON_LEFT) and not (self.place):
            #Run selecting object code here
            self.selected_object = arcade.get_sprites_at_point((x, y), self.objects)

            # Check to see if we actually clicked on an object: If yes make that the current object, if not, make None
            if self.selected_object:
                self.selected_object = self.selected_object[-1] #-1 to grab the top-most object
            else:
                self.selected_object = None
            print(self.selected_object)



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

        if (symbol == arcade.key.I) and (self.selected_object is not None):
            print(self.selected_object.pos)

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


        # Deletes selected object
        if (symbol == arcade.key.DELETE) and (self.selected_object is not None):
            self.selected_object.kill()
            self.selected_object = None


        #Show Vectors
        if symbol == arcade.key.V:
            self.draw_vecs = not self.draw_vecs


        #Spawn random earthlikes
        if symbol == arcade.key.M:
            for _ in range(5):
                self.objects.append(Entity(25_000_000_000_000_000, rand.randint(0, self.width), rand.randint(0, self.height), "media/earthlike.png"))


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
