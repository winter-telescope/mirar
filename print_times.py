from astroplan import observer
from astroplan import Observer
from datetime import datetime
from astropy.time import Time

palomar = Observer.at_site('Palomar')
t = Time(datetime.utcnow())

twilight_6 = palomar.twilight_evening_civil(t,which='next')
twilight_12 = palomar.twilight_evening_nautical(t,which='next')
twilight_18 = palomar.twilight_evening_astronomical(t,which='next')
sunset = palomar.sun_set_time(t,which='next')

print('Sunset at ',sunset.isot)
print('6-degrees at ',twilight_6.isot)
print('12-degrees at ',twilight_12.isot)
print('18-degrees twilight at ',twilight_18.isot)