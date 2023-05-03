import os
from string import Template

WRAPPER_NAMES = ["wrapper.py", "wrapper.r"]
ENVIRONMENT_YAML_FILE = "environment.yaml"


TEMPLATE_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(TEMPLATE_DIR)

OUT_DIR = os.path.join(BASE_DIR, "_docs")
WRAPPERS_DIR = os.path.join(BASE_DIR, "wrappers")
DOCS_DIR = os.path.join(OUT_DIR, "docs")
MKDOCS_YML = os.path.join(OUT_DIR, "mkdocs.yml")

TEMPLATE_MKDOCS = os.path.join(TEMPLATE_DIR, "mkdocs.template")
TEMPLATE_WRAPPER = os.path.join(TEMPLATE_DIR, "wrapper.md")


def generate_documentation():
    os.makedirs(DOCS_DIR, exist_ok=True)

    navigation_list: list[str] = []
    for tool, tool_dir in list_dirs_in_dir(WRAPPERS_DIR):
        navigation_list.append(first_level(tool))

        for command, path in list_dirs_in_dir(tool_dir):
            md_file = f"{tool}_{command}.md"
            navigation_list.append(second_level(command, md_file))

            md_code = markdown_for_command(tool, command, path)

            md_fullpath = os.path.join(DOCS_DIR, md_file)
            write_to_file(md_code, md_fullpath)

    mkdocs_yml_code = populate_mkdocs_with_navigation(navigation_list)
    write_to_file(mkdocs_yml_code, MKDOCS_YML)


def list_dirs_in_dir(listing_dir: str) -> list[tuple[str, str]]:
    return [
        (name, full_path)
        for name in os.listdir(listing_dir)
        if os.path.isdir(full_path := os.path.join(listing_dir, name))
        and not name.startswith((".", "_"))
    ]


def first_level(name: str) -> str:
    return f"  - {name}:"


def second_level(name: str, ref_file: str) -> str:
    return f"    - {name}: {ref_file}"


def write_to_file(content: str, path: str):
    with open(path, "w") as f:
        f.write(content)


def populate_mkdocs_with_navigation(navigation_list: list[str]):
    with open(TEMPLATE_MKDOCS, "r") as template_mkdocs:
        mkdocs_yml_code = Template(template_mkdocs.read()).substitute(
            LIST="\n".join(navigation_list),
        )
    return mkdocs_yml_code


def find_wrapper_code(path: str) -> str:
    for wrapper in WRAPPER_NAMES:
        possible_wrapper = os.path.join(path, wrapper)
        if os.path.exists(possible_wrapper):
            return possible_wrapper
    raise ValueError(f'Could not find wrapper for "{path}"')


def find_command_code_language(path: str):
    extension = os.path.splitext(path)
    if extension[1] == ".py":
        return "python"
    elif extension[1] == ".r":
        return "r"
    else:
        raise ValueError(f'Could not determine command code language for "{path}"')


def markdown_for_command(tool: str, command: str, path: str) -> str:
    wrapper_code_path = find_wrapper_code(path)
    wrapper_code = read_file(wrapper_code_path)

    env_yaml_path = os.path.join(path, ENVIRONMENT_YAML_FILE)
    env_yaml = read_file(env_yaml_path)

    snakefile_path = os.path.join(path, "test", "Snakefile")
    snakefile = read_file(snakefile_path)

    readme_path = os.path.join(path, "README.md")
    readme_content = read_file(readme_path)

    with open(TEMPLATE_WRAPPER, "r") as md_template:
        md_code = Template(md_template.read()).substitute(
            command_description=readme_content,
            command_code=wrapper_code,
            command_language=find_command_code_language(wrapper_code_path),
            environment_yaml=env_yaml,
            snakefile_code=snakefile,
        )
    return md_code


def read_file(file: str) -> str:
    with open(file, "r") as f:
        file_content = f.read().rstrip()
    return file_content


if __name__ == "__main__":
    generate_documentation()
