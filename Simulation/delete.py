import pandas as pd
import requests
import json
import cloudscraper

url = "https://api.cryptonator.com/api/full/btc-USD"

if __name__ == "__main__":
    scraper = cloudscraper.create_scraper()
    # json file moet geskryf word na current repository en gelees word
    with open('data.json', 'w') as f:
        json.dump(scraper.get(url).json(), f)
    
    df = pd.read_json('data.json')
    print(df)
    markets = df['ticker']['markets']

    mark_df = pd.DataFrame(columns=['market', 'price', 'volume'], index = [0])

    market_list = []

    for market in markets:
        for col in mark_df.columns:
            mark_df[col][0] = market[col]
        market_list.append(mark_df.copy())

    markets = pd.concat(market_list)
    markets.reset_index(drop=True, inplace=True)

    print(markets)