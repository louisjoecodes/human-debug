import React, { useEffect, useRef } from 'react';
import EditorJS from '@editorjs/editorjs';
import Header from '@editorjs/header';
import List from '@editorjs/list';
import Paragraph from '@editorjs/paragraph';

interface EditorProps {
    defaultValue?: any;
    onChange?: (data: any) => void;
}

const Editor: React.FC<EditorProps> = ({ defaultValue, onChange }) => {
    const editorRef = useRef<EditorJS | null>(null);

    useEffect(() => {
        if (!editorRef.current) {
            const editor = new EditorJS({
                holder: 'editorjs',
                tools: {
                    header: Header,
                    list: List,
                    paragraph: Paragraph,
                },
                data: defaultValue,
                onChange: async () => {
                    const content = await editorRef.current?.save();
                    onChange && onChange(content);
                },
            });

            editorRef.current = editor;
        }

        return () => {
            if (editorRef.current && editorRef.current.destroy) {
                editorRef.current.destroy();
            }
        };
    }, []);

    return <div id="editorjs" className="prose max-w-none" />;
};

export default Editor;