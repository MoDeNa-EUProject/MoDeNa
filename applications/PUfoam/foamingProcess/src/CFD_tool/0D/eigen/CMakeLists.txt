file(GLOB eigen_SRC
    "src/*.f"
)

add_library(eigen ${eigen_SRC})

