import os
import folium
import pandas as pd
import json
from folium import GeoJson
from folium import IFrame

haiti = json.load(open('haiti.geojson', 'r+', encoding='utf-8'))
water = json.load(open('hti_dom_watcrsl_rvr_osm.geojson', 'r+', encoding='utf=8'))

attr = ('&copy; <a href="http://www.openstreetmap.org/copyright">OpenStreetMap</a> '
        'contributors, &copy; <a href="http://cartodb.com/attributions">CartoDB</a>')
tiles = 'http://{s}.basemaps.cartocdn.com/light_nolabels/{z}/{x}/{y}.png'

m = folium.Map(location=[19, -72], zoom_start=9, attr=attr, tiles=tiles)

folium.GeoJson(haiti, style_function=lambda x:
	{'color' : x['properties']['stroke'], 'weight' :
	x['properties']['stroke-width'], 'strokeOpacity': x['properties']['stroke-opacity'],
	'fillColor' : x['properties']['fill'], 'fillOpacity': x['properties']['fill-opacity']}
	).add_to(m)

folium.GeoJson(water, style_function=lambda x:
  {'color' : '#58bbff', 'stroke' : '#58bbff', 'fill' : '#58bbff', 'strokeWidth' : '0.25', 
  'strokeOpacity': '0.75', 'fillColor' : '#58bbff', 'fillOpacity': '0.25'}
  ).add_to(m)


m.save('index.html')