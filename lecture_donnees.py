#-*- coding: utf-8 -*-
from math import sqrt
import numpy as np
from structures import Point, Velocity

__author__ = "spraakforskaren"


# class Mesures_balise:
#     def __init__(self):



def read_acquired_data_from_beacons(nom_fichier):
    with open(nom_fichier, "r") as f:
        raw_data = f.read()
    return [line.replace("#80#", "").split(";") for line in raw_data.split("\n") if line != ""]


def change_unity(n):
    """

    :param n: is a string reprenting a number
    :return:
    """
    return str(float(n)/1000000*340000)


def permute_measure(line, permutation):
    # return [change_unity(line[permutation[0]]), change_unity(line[permutation[1]]), change_unity(line[permutation[2]])]
    t1 = change_unity(line[permutation[0]])
    t2 = change_unity(line[permutation[1]])
    t3 = change_unity(line[permutation[2]])
    return [str(float(t1)-float(t2)), str(float(t1)-float(t3)), str(float(t2)-float(t3))]


def convert_system_for(data, nom_fichier, name, permutation):
    with open("converted_"+name+"_"+nom_fichier, "w") as f:
        for line in data:
            f.write(";".join(permute_measure(line, permutation))+"\n")


def convert_to_clement_system(data, nom_fichier):
    permutation = [0, 2, 1]
    name = "clement"
    convert_system_for(data, nom_fichier, name, permutation)


def convert_to_ivan_system(data, nom_fichier):
    permutation = [1, 2, 0]
    name = "ivan"
    convert_system_for(data, nom_fichier, name, permutation)


def script_convert_to_clement_system(nom_fichier):
    data = read_acquired_data_from_beacons(nom_fichier)
    convert_to_clement_system(data, nom_fichier)


def script_convert_to_ivan_system(nom_fichier):
    data = read_acquired_data_from_beacons(nom_fichier)
    convert_to_ivan_system(data, nom_fichier)


if __name__ == "__main__":
    script_convert_to_clement_system("acquisition_1.txt")
    script_convert_to_clement_system("acquisition_2.txt")
    script_convert_to_ivan_system("acquisition_1.txt")
    script_convert_to_ivan_system("acquisition_2.txt")
