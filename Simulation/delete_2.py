from luno_python.client import Client

if __name__ == '__main__':

    id = "buye79xc5dtdw"

    password = "mtuPk8XNQ7RVIUnClv5PnOSIwTadEfxiy4cCXKH2QK0"

    c = Client(api_key_id=id, api_key_secret=password)

    print(c.get_balances())