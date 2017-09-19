from TrusightOne.jsonhandler import JsonHandler

# import pickle

json_dir = "/media/kasutaja/data/NGS_data/panels_json_main/"


def main():
    handler = JsonHandler(json_dir)
    panels = handler.get_all_panels(external=False)
    for panel in panels:
        print(panel)


if __name__ == "__main__":
    main()
