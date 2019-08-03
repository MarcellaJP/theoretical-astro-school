set_import_module_path( strcat(".", char(path_get_delimiter()), get_import_module_path()) );
import("simplejet");
exit(0);
