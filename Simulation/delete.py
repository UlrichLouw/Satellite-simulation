import pandas as pd
import requests
import json
import cloudscraper

url = "https://api.cryptonator.com/api/full/btc-usd"

if __name__ == "__main__":
    scraper = cloudscraper.create_scraper()
    # json file moet geskryf word na current repository en gelees word
    df = pd.from_json(scraper.get(url).json())
    print(df)