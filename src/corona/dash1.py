"""Dashboards"""
##
from corona.utils import *
from corona.maths import *
from corona.model import *
from corona.plotting import *

import mpl_tools

## Download data
import requests
import requests_cache
# Cache persitence/expiration
from datetime import timedelta
expire_after = timedelta(hours=1)

def get_json(URL):
    # from urllib.parse import urlparse
    # shortened = urlparse(URL).netloc
    import pathlib
    shortened = pathlib.Path(URL).name
    requests_cache.install_cache(shortened, expire_after=expire_after)
    data = requests.get(URL)
    source = "cache" if data.from_cache else "internet"
    print(f"{shortened} from {source}.")
    return data.json()

def _c(c):
    if c=="United States": return "US"
    if c=="United Kingdom": return "UK"
    return c

covid19 = get_json('https://pomber.github.io/covid19/timeseries.json')
covid19 = {_c(k):v for k,v in covid19.items()}

# Population data
pops = get_json('https://raw.githubusercontent.com/samayo/country-json/master/src/country-by-population.json')
byCountry = {_c(d['country']): {"nPop":int(d['population'] or -1)} for d in pops}

# As dict
byDate = {dct["date"]: {k:v for k,v in dct.items() if k!="date"}  for dct in covid19["Norway"]}
# IPython.lib.pretty.pretty prints better than default, but this is still better:
byDate = JsonDict(byDate)

# As dataframes
import pandas as pd
def get_df(country):
    d  = covid19[country]
    df = pd.DataFrame.from_records(d)
    df["date"] = pd.to_datetime(df["date"])
    return df.set_index("date")
dfs = {country: get_df(country) for country in covid19}

# Examples
# df = dfs["Norway"]
# column = df["deaths"]
# df[(1<=df["deaths"]) & (df["deaths"]<=3)]
# df[0!=df["deaths"]]
# Accessing rows:
# print(dfs["Norway"][-5:])
# print(dfs["Norway"].tail())
# print(dfs["Norway"].loc["2020-3-29":])
# print(dfs["Norway"]["2020-3-29":])


## Plot all series for 1 country
c = "Norway"
fig, ax = mpl_tools.freshfig(1)
dfs[c].plot(ax=ax)
mpl_tools.add_log_toggler(ax)
ax.set_title(c)

## Plot 1 serie for multiple countries
fig, ax = mpl_tools.freshfig(2)
series = "deaths"
countries = "Norway", "Sweden", "France", "UK", "US"
for c in countries:
    dfs[c][series].plot(ax=ax,label=c)
ax.legend()
ax.set_title(series)
mpl_tools.add_log_toggler(ax)


## 
fig, ax = mpl_tools.freshfig(3)
deaths_intercept = 1 / 10**6 # death (density) for which lines should intersect
for c in countries:
    # Compute death density
    nPop = byCountry[c]["nPop"]
    deaths = dfs[c]["deaths"].to_numpy()/nPop
    # Compute date offset
    day_range = arange(len(deaths))
    nDay = np.interp(deaths_intercept, deaths, day_range)
    # Plot with offset
    ax.plot(day_range - nDay, deaths, label=c)
    # Store offset
    byCountry[c]["lag"] = nDay
ax.plot([0,0],ax.get_ylim(), "k--", lw=1,label="_nolegend_")
ax.set_ylabel("Deaths per capita (cDeath)")
ax.set_title(f"Day 0 defined by cDeath(t)={deaths_intercept}")
ax.set_xlabel(f"Days (t)")
ax.legend()
mpl_tools.add_log_toggler(ax)
##

        


##
